#include <iostream>
#include <llvm/ADT/SmallVector.h>
#include <llvm/Analysis/CGSCCPassManager.h>
#include <llvm/Analysis/LoopAnalysisManager.h>
#include <llvm/Bitcode/BitcodeWriter.h>
#include <llvm/Passes/PassBuilder.h>
#include <llvm/Support/CommandLine.h>
#include <llvm/Support/Error.h>
#include <llvm/Support/InitLLVM.h>
#include <llvm/Support/MemoryBuffer.h>
#include <llvm/Support/TargetSelect.h>
#include <llvm/Transforms/Utils.h>
#include <marco/ast/Parser.h>
#include <marco/ast/Passes.h>
#include <marco/codegen/CodeGen.h>
#include <marco/codegen/Passes.h>
#include <marco/utils/VariableFilter.h>
#include <mlir/Conversion/Passes.h>
#include <mlir/Dialect/LLVMIR/LLVMDialect.h>
#include <mlir/Dialect/StandardOps/IR/Ops.h>
#include <mlir/Dialect/StandardOps/Transforms/Passes.h>
#include <mlir/ExecutionEngine/OptUtils.h>
#include <mlir/IR/BuiltinOps.h>
#include <mlir/IR/MLIRContext.h>
#include <mlir/Pass/PassManager.h>
#include <mlir/Target/LLVMIR/Dialect/LLVMIR/LLVMToLLVMIRTranslation.h>
#include <mlir/Target/LLVMIR/Dialect/OpenMP/OpenMPToLLVMIRTranslation.h>
#include <mlir/Target/LLVMIR/Export.h>
#include <mlir/Transforms/Passes.h>

#include "clang/Driver/Driver.h"
#include "clang/Basic/Diagnostic.h"
#include "clang/Basic/DiagnosticIDs.h"
#include "clang/Basic/DiagnosticOptions.h"
#include "clang/Driver/Compilation.h"
#include "llvm/ADT/ArrayRef.h"
#include "llvm/ADT/IntrusiveRefCntPtr.h"
#include "llvm/Option/ArgList.h"
#include "llvm/Support/Host.h"
#include "llvm/Support/InitLLVM.h"
#include "llvm/Support/VirtualFileSystem.h"

#include <marco/frontend/CompilerInstance.h>
#include <marco/frontend/CompilerInvocation.h>
#include <marco/frontend/TextDiagnosticBuffer.h>
#include <marco/frontendTool/Utils.h>
#include <clang/Driver/DriverDiagnostic.h>
#include <llvm/Option/Arg.h>
#include <llvm/Option/ArgList.h>
#include <llvm/Option/OptTable.h>

#include <marco/frontend/CompilerInvocation.h>
#include <marco/frontend/TextDiagnosticPrinter.h>

using namespace llvm;
using namespace marco;
using namespace marco::frontend;
using namespace std;

extern int mc1_main(llvm::ArrayRef<const char *> argv, const char *argv0);
extern int cc1_main(ArrayRef<const char*> argv, const char* argv0, void* mainAddr);

static std::string getExecutablePath(const char* argv0)
{
  // This just needs to be some symbol in the binary
  void* p = (void*) (intptr_t) getExecutablePath;
  return llvm::sys::fs::getMainExecutable(argv0, p);
}

static int executeMC1Tool(llvm::SmallVectorImpl<const char *>& argv)
{
  llvm::StringRef tool = argv[1];

  if (tool == "-mc1") {
    return mc1_main(makeArrayRef(argv).slice(2), argv[0]);
  }

  // Reject unknown tools
  llvm::errs() << "error: unknown integrated tool '" << tool << "'. "
               << "Valid tools include '-mc1'.\n";

  return 1;
}

static int executeCC1Tool(SmallVectorImpl<const char*>& argv) {
  // If we call the cc1 tool from the clangDriver library (through
  // Driver::CC1Main), we need to clean up the options usage count. The options
  // are currently global, and they might have been used previously by the
  // driver.
  llvm::cl::ResetAllOptionOccurrences();

  llvm::BumpPtrAllocator A;
  llvm::StringSaver Saver(A);
  llvm::cl::ExpandResponseFiles(Saver, &llvm::cl::TokenizeGNUCommandLine, argv, false);

  StringRef tool = argv[1];

  void* GetExecutablePathVP = (void*)(intptr_t)getExecutablePath;

  if (tool == "-cc1") {
    return cc1_main(makeArrayRef(argv).slice(1), argv[0], GetExecutablePathVP);
  }

  // Reject unknown tools
  llvm::errs() << "error: unknown integrated tool '" << tool << "'. "
               << "Valid tools include '-cc1'.\n";

  return 1;
}

// This lets us create the DiagnosticsEngine with a properly-filled-out
// DiagnosticOptions instance
static clang::DiagnosticOptions *CreateAndPopulateDiagOpts(
    llvm::ArrayRef<const char *> argv) {
  auto *diagOpts = new clang::DiagnosticOptions;

  // Ignore missingArgCount and the return value of ParseDiagnosticArgs.
  // Any errors that would be diagnosed here will also be diagnosed later,
  // when the DiagnosticsEngine actually exists.
  unsigned missingArgIndex, missingArgCount;

  llvm::opt::InputArgList args = clang::driver::getDriverOptTable().ParseArgs(
      argv.slice(1), missingArgIndex, missingArgCount);

  return diagOpts;
}

int main(int argc, const char** argv)
{
  // Initialize variables to call the driver
  llvm::InitLLVM x(argc, argv);

  llvm::SmallVector<const char *, 256> args(argv, argv + argc);

  // Check if MARCO is in the frontend mode
  auto firstArg = std::find_if(args.begin() + 1, args.end(), [](const char* arg) {
    return arg != nullptr;
  });

  if (firstArg != args.end()) {
    if (llvm::StringRef(*firstArg).startswith("-mc1")) {
      return executeMC1Tool(args);
    }

    if (llvm::StringRef(*firstArg).startswith("-cc1")) {
      return executeCC1Tool(args);
    }
  }

  // Not in the frontend mode. Continue in the compiler driver mode.
  llvm::InitializeAllTargets();
  llvm::InitializeAllTargetMCs();
  llvm::InitializeAllAsmPrinters();

  // Create the diagnostics engine for the driver
  llvm::IntrusiveRefCntPtr<clang::DiagnosticOptions> diagOpts = CreateAndPopulateDiagOpts(args);
  llvm::IntrusiveRefCntPtr<clang::DiagnosticIDs> diagID(new clang::DiagnosticIDs());

  marco::frontend::TextDiagnosticPrinter* diagClient =
      new marco::frontend::TextDiagnosticPrinter(llvm::errs(), &*diagOpts);

  diagClient->set_prefix(std::string(llvm::sys::path::stem(getExecutablePath(args[0]))));

  clang::DiagnosticsEngine diags(diagID, &*diagOpts, diagClient);

  // Prepare the driver
  auto triple = llvm::sys::getDefaultTargetTriple();
  clang::driver::ParsedClangName targetAndMode(triple, "cc1", "--driver-mode=cc1", true);
  std::string driverPath = getExecutablePath(args[0]);

  clang::driver::Driver theDriver(driverPath, triple, diags, "MARCO");
  theDriver.setTargetAndMode(targetAndMode);
  std::unique_ptr<clang::driver::Compilation> c(theDriver.BuildCompilation(args));
  llvm::SmallVector<std::pair<int, const clang::driver::Command*>, 4> failingCommands;

  // Run the driver
  int res = 1;
  bool isCrash = false;
  res = theDriver.ExecuteCompilation(*c, failingCommands);

  for (const auto& p: failingCommands) {
    int CommandRes = p.first;
    const clang::driver::Command* failingCommand = p.second;
    if (!res) {
      res = CommandRes;
    }

    // If result status is < 0 (e.g. when sys::ExecuteAndWait returns -1),
    // then the driver command signalled an error. On Windows, abort will
    // return an exit code of 3. In these cases, generate additional diagnostic
    // information if possible.
    isCrash = CommandRes < 0;
    #ifdef _WIN32
    isCrash |= CommandRes == 3;
    #endif

    if (isCrash) {
      theDriver.generateCompilationDiagnostics(*c, *failingCommand);
      break;
    }
  }

  diags.getClient()->finish();

  // If we have multiple failing commands, we return the result of the first
  // failing command.
  return res;
}
