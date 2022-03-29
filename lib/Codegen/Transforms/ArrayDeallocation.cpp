#include "marco/Codegen/Transforms/ArrayDeallocation.h"
#include "marco/Dialect/Modelica/ModelicaDialect.h"
#include "marco/Codegen/Transforms/Model/ModelConverter.h"
#include "mlir/Transforms/BufferUtils.h"
#include "mlir/Transforms/Passes.h"

using namespace ::marco::codegen;
using namespace ::mlir::modelica;

namespace
{
  class ArrayDeallocation : public mlir::BufferPlacementTransformationBase
  {
    public:
      ArrayDeallocation(mlir::Operation* op)
          : BufferPlacementTransformationBase(op),
            postDominators(op)
      {
      }

      void deallocate() const
      {
        for (const mlir::BufferPlacementAllocs::AllocEntry& entry : allocs) {
          mlir::Value alloc = std::get<0>(entry);

          if (auto arrayType = alloc.getType().dyn_cast<ArrayType>(); arrayType && !arrayType.hasConstantShape()) {
            bool isStored = llvm::any_of(alloc.getUsers(), [&](const auto& op) {
              if (auto memberStoreOp = mlir::dyn_cast<MemberStoreOp>(op)) {
                return memberStoreOp.value() == alloc;
              }

              return false;
            });

            if (isStored) {
              continue;
            }
          }

          auto aliasesSet = aliases.resolve(alloc);
          assert(!aliasesSet.empty() && "Must contain at least one alias");

          // Determine the actual block to place the dealloc and get liveness information.
          mlir::Block* placementBlock = findCommonDominator(alloc, aliasesSet, postDominators);
          const mlir::LivenessBlockInfo* livenessInfo = liveness.getLiveness(placementBlock);

          // We have to ensure that the dealloc will be after the last use of all
          // aliases of the given value. We first assume that there are no uses in
          // the placementBlock and that we can safely place the dealloc at the
          // beginning.
          mlir::Operation* endOperation = &placementBlock->front();

          // Iterate over all aliases and ensure that the endOperation will point
          // to the last operation of all potential aliases in the placementBlock.
          for (mlir::Value alias : aliasesSet) {
            // Ensure that the start operation is at least the defining operation of the current alias to avoid invalid placement of deallocs for aliases without any uses.
            mlir::Operation* beforeOp = endOperation;

            if (alias.getDefiningOp() && !(beforeOp = placementBlock->findAncestorOpInBlock(*alias.getDefiningOp()))) {
              continue;
            }

            mlir::Operation* aliasEndOperation = livenessInfo->getEndOperation(alias, beforeOp);

            // Check whether the aliasEndOperation lies in the desired block and
            // whether it is behind the current endOperation. If yes, this will be
            // the new endOperation.

            if (aliasEndOperation->getBlock() == placementBlock && endOperation->isBeforeInBlock(aliasEndOperation)) {
              endOperation = aliasEndOperation;
            }
          }

          // endOperation is the last operation behind which we can safely store
          // the dealloc taking all potential aliases into account.

          // If there is an existing dealloc, move it to the right place.
          mlir::Operation* deallocOperation = std::get<1>(entry);

          if (deallocOperation) {
            deallocOperation->moveAfter(endOperation);
          } else {
            // If the Dealloc position is at the terminator operation of the
            // block, then the value should escape from a deallocation.

            mlir::Operation* nextOp = endOperation->getNextNode();

            if (!nextOp) {
              continue;
            }

            // If there is no dealloc node, insert one in the right place.
            mlir::OpBuilder builder(nextOp);
            builder.create<FreeOp>(alloc.getLoc(), alloc);
          }
        }
      }

    private:
      mlir::PostDominanceInfo postDominators;
  };

  class ArrayDeallocationPass : public mlir::PassWrapper<ArrayDeallocationPass, mlir::OperationPass<mlir::ModuleOp>>
  {
    public:
      void runOnOperation() override
      {
        getOperation().walk([](mlir::FuncOp op) {
          if (op.getName() != ModelConverter::initFunctionName) {
            ArrayDeallocation deallocation(op);
            deallocation.deallocate();
          }
        });

        getOperation().walk([](FunctionOp op) {
          ArrayDeallocation deallocation(op);
          deallocation.deallocate();
        });
      }
  };
}

namespace marco::codegen
{
  std::unique_ptr<mlir::Pass>createArrayDeallocationPass()
  {
    return std::make_unique<ArrayDeallocationPass>();
  }
}
