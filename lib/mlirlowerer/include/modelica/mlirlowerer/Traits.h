#pragma once

#include <mlir/IR/OpDefinition.h>
#include <mlir/IR/BlockAndValueMapping.h>

namespace detail
{
	struct EquationInterfaceTraits
	{
		struct Concept
		{
			Concept() = default;
			Concept(const Concept& other) = default;
			Concept(Concept&& other) = default;
			Concept& operator=(Concept&& other) = default;
			virtual ~Concept() = default;
			Concept& operator=(const Concept& other) = default;

			virtual mlir::Block* body(mlir::Operation* op) const = 0;
			virtual mlir::ValueRange inductions(mlir::Operation* op) const = 0;
			virtual mlir::Value induction(mlir::Operation* op, size_t index) const = 0;
			virtual long inductionIndex(mlir::Operation* op, mlir::Value induction) const = 0;
			virtual mlir::ValueRange lhs(mlir::Operation* op) const = 0;
			virtual mlir::ValueRange rhs(mlir::Operation* op) const = 0;
		};

		template<typename ConcreteOp>
		struct Model : public Concept
		{
			mlir::Block* body(mlir::Operation* op) const final
			{
				return mlir::cast<ConcreteOp>(op).body();
			}

			mlir::ValueRange inductions(mlir::Operation* op) const final
			{
				return mlir::cast<ConcreteOp>(op).inductions();
			}

			mlir::Value induction(mlir::Operation* op, size_t index) const final
			{
				return mlir::cast<ConcreteOp>(op).induction(index);
			}

			long inductionIndex(mlir::Operation* op, mlir::Value induction) const final
			{
				return mlir::cast<ConcreteOp>(op).inductionIndex(induction);
			}

			mlir::ValueRange lhs(mlir::Operation* op) const final
			{
				return mlir::cast<ConcreteOp>(op).lhs();
			}

			mlir::ValueRange rhs(mlir::Operation* op) const final
			{
				return mlir::cast<ConcreteOp>(op).rhs();
			}
		};

		template<typename ConcreteOp>
		class FallbackModel : public Concept
		{
			public:
			FallbackModel() = default;

			mlir::Block* body(mlir::Operation* op) const final
			{
				return mlir::cast<ConcreteOp>(op).body();
			}

			mlir::ValueRange inductions(mlir::Operation* op) const final
			{
				return mlir::cast<ConcreteOp>(op).inductions();
			}

			mlir::Value induction(mlir::Operation* op, size_t index) const final
			{
				return mlir::cast<ConcreteOp>(op).induction(index);
			}

			long inductionIndex(mlir::Operation* op, mlir::Value induction) const final
			{
				return mlir::cast<ConcreteOp>(op).inductionIndex(induction);
			}

			mlir::ValueRange lhs(mlir::Operation* op) const final
			{
				return mlir::cast<ConcreteOp>(op).lhs();
			}

			mlir::ValueRange rhs(mlir::Operation* op) const final
			{
				return mlir::cast<ConcreteOp>(op).rhs();
			}
		};
	};
}

class EquationInterface : public mlir::OpInterface<EquationInterface, detail::EquationInterfaceTraits>
{
	public:
	using mlir::OpInterface<EquationInterface, detail::EquationInterfaceTraits>::OpInterface;

	mlir::Block* body()
	{
		return getImpl()->body(getOperation());
	}

	mlir::ValueRange inductions()
	{
		return getImpl()->inductions(getOperation());
	}

	mlir::Value induction(size_t index)
	{
		return getImpl()->induction(getOperation(), index);
	}

	long inductionIndex(mlir::Value induction)
	{
		return getImpl()->inductionIndex(getOperation(), induction);
	}

	mlir::ValueRange lhs()
	{
		return getImpl()->lhs(getOperation());
	}

	mlir::ValueRange rhs()
	{
		return getImpl()->rhs(getOperation());
	}
};

namespace detail
{
	struct InvertibleInterfaceTraits
	{
		struct Concept
		{
			Concept() = default;
			Concept(const Concept& other) = default;
			Concept(Concept&& other) = default;
			Concept& operator=(Concept&& other) = default;
			virtual ~Concept() = default;
			Concept& operator=(const Concept& other) = default;

			// TODO: keep ValueRange or switch to Value?
			virtual mlir::LogicalResult invert(mlir::Operation* op, mlir::OpBuilder& builder, unsigned int argumentIndex, mlir::ValueRange currentResult) const = 0;
		};

		template <typename ConcreteOp>
		struct Model : public Concept
		{
			mlir::LogicalResult invert(mlir::Operation* op, mlir::OpBuilder& builder, unsigned int argumentIndex, mlir::ValueRange currentResult) const final
			{
				return mlir::cast<ConcreteOp>(op).invert(builder, argumentIndex, currentResult);
			}
		};

		template<typename ConcreteOp>
		class FallbackModel : public Concept
		{
			public:
			FallbackModel() = default;

			mlir::LogicalResult invert(mlir::Operation* op, mlir::OpBuilder& builder, unsigned int argumentIndex, mlir::ValueRange currentResult) const final
			{
				return mlir::cast<ConcreteOp>(op).invert(builder, argumentIndex, currentResult);
			}
		};
	};
}

class InvertibleInterface : public mlir::OpInterface<InvertibleInterface, detail::InvertibleInterfaceTraits>
{
	public:
	using OpInterface<InvertibleInterface, detail::InvertibleInterfaceTraits>::OpInterface;

	mlir::LogicalResult invert(mlir::OpBuilder& builder, unsigned int argumentIndex, mlir::ValueRange currentResult)
	{
		return getImpl()->invert(getOperation(), builder, argumentIndex, currentResult);
	}
};

namespace detail
{
	struct DistributableInterfaceTraits
	{
		struct Concept
		{
			Concept() = default;
			Concept(const Concept& other) = default;
			Concept(Concept&& other) = default;
			Concept& operator=(Concept&& other) = default;
			virtual ~Concept() = default;
			Concept& operator=(const Concept& other) = default;

			virtual mlir::Value distribute(mlir::Operation* op, mlir::OpBuilder& builder) const = 0;
		};

		template <typename ConcreteOp>
		struct Model : public Concept
		{
			mlir::Value distribute(mlir::Operation* op, mlir::OpBuilder& builder) const final
			{
				return mlir::cast<ConcreteOp>(op).distribute(builder);
			}
		};

		template<typename ConcreteOp>
		class FallbackModel : public Concept
		{
			public:
			FallbackModel() = default;

			mlir::Value distribute(mlir::Operation* op, mlir::OpBuilder& builder) const final
			{
				return mlir::cast<ConcreteOp>(op).distribute(builder);
			}
		};
	};
}

class DistributableInterface : public mlir::OpInterface<DistributableInterface, detail::DistributableInterfaceTraits>
{
	public:
	using OpInterface<DistributableInterface, detail::DistributableInterfaceTraits>::OpInterface;

	mlir::Value distribute(mlir::OpBuilder& builder)
	{
		return getImpl()->distribute(getOperation(), builder);
	}
};

namespace detail
{
	struct NegateOpDistributionInterfaceTraits
	{
		struct Concept
		{
			Concept() = default;
			Concept(const Concept& other) = default;
			Concept(Concept&& other) = default;
			Concept& operator=(Concept&& other) = default;
			virtual ~Concept() = default;
			Concept& operator=(const Concept& other) = default;

			virtual mlir::Value distributeNegateOp(mlir::Operation* op, mlir::OpBuilder& builder, mlir::Type resultType) const = 0;
		};

		template <typename ConcreteOp>
		struct Model : public Concept
		{
			mlir::Value distributeNegateOp(mlir::Operation* op, mlir::OpBuilder& builder, mlir::Type resultType) const final
			{
				return mlir::cast<ConcreteOp>(op).distributeNegateOp(builder, resultType);
			}
		};

		template<typename ConcreteOp>
		class FallbackModel : public Concept
		{
			public:
			FallbackModel() = default;

			mlir::Value distributeNegateOp(mlir::Operation* op, mlir::OpBuilder& builder, mlir::Type resultType) const final
			{
				return mlir::cast<ConcreteOp>(op).distributeNegateOp(builder, resultType);
			}
		};
	};
}

class NegateOpDistributionInterface : public mlir::OpInterface<NegateOpDistributionInterface, detail::NegateOpDistributionInterfaceTraits>
{
	public:
	using OpInterface<NegateOpDistributionInterface, detail::NegateOpDistributionInterfaceTraits>::OpInterface;

	mlir::Value distributeNegateOp(mlir::OpBuilder& builder, mlir::Type resultType)
	{
		return getImpl()->distributeNegateOp(getOperation(), builder, resultType);
	}
};

namespace detail
{
	struct MulOpDistributionInterfaceTraits
	{
		struct Concept
		{
			Concept() = default;
			Concept(const Concept& other) = default;
			Concept(Concept&& other) = default;
			Concept& operator=(Concept&& other) = default;
			virtual ~Concept() = default;
			Concept& operator=(const Concept& other) = default;

			virtual mlir::Value distributeMulOp(mlir::Operation* op, mlir::OpBuilder& builder, mlir::Type resultType, mlir::Value value) const = 0;
		};

		template <typename ConcreteOp>
		struct Model : public Concept
		{
			mlir::Value distributeMulOp(mlir::Operation* op, mlir::OpBuilder& builder, mlir::Type resultType, mlir::Value value) const final
			{
				return mlir::cast<ConcreteOp>(op).distributeMulOp(builder, resultType, value);
			}
		};

		template<typename ConcreteOp>
		class FallbackModel : public Concept
		{
			public:
			FallbackModel() = default;

			mlir::Value distributeMulOp(mlir::Operation* op, mlir::OpBuilder& builder, mlir::Type resultType, mlir::Value value) const final
			{
				return mlir::cast<ConcreteOp>(op).distributeMulOp(builder, resultType, value);
			}
		};
	};
}

class MulOpDistributionInterface : public mlir::OpInterface<MulOpDistributionInterface, detail::MulOpDistributionInterfaceTraits>
{
	public:
	using OpInterface<MulOpDistributionInterface, detail::MulOpDistributionInterfaceTraits>::OpInterface;

	mlir::Value distributeMulOp(mlir::OpBuilder& builder, mlir::Type resultType, mlir::Value value)
	{
		return getImpl()->distributeMulOp(getOperation(), builder, resultType, value);
	}
};

namespace detail
{
	struct DivOpDistributionInterfaceTraits
	{
		struct Concept
		{
			Concept() = default;
			Concept(const Concept& other) = default;
			Concept(Concept&& other) = default;
			Concept& operator=(Concept&& other) = default;
			virtual ~Concept() = default;
			Concept& operator=(const Concept& other) = default;

			virtual mlir::Value distributeDivOp(mlir::Operation* op, mlir::OpBuilder& builder, mlir::Type resultType, mlir::Value value) const = 0;
		};

		template <typename ConcreteOp>
		struct Model : public Concept
		{
			mlir::Value distributeDivOp(mlir::Operation* op, mlir::OpBuilder& builder, mlir::Type resultType, mlir::Value value) const final
			{
				return mlir::cast<ConcreteOp>(op).distributeDivOp(builder, resultType, value);
			}
		};

		template<typename ConcreteOp>
		class FallbackModel : public Concept
		{
			public:
			FallbackModel() = default;

			mlir::Value distributeDivOp(mlir::Operation* op, mlir::OpBuilder& builder, mlir::Type resultType, mlir::Value value) const final
			{
				return mlir::cast<ConcreteOp>(op).distributeDivOp(builder, resultType, value);
			}
		};
	};
}

class DivOpDistributionInterface : public mlir::OpInterface<DivOpDistributionInterface, detail::DivOpDistributionInterfaceTraits>
{
	public:
	using OpInterface<DivOpDistributionInterface, detail::DivOpDistributionInterfaceTraits>::OpInterface;

	mlir::Value distributeDivOp(mlir::OpBuilder& builder, mlir::Type resultType, mlir::Value value)
	{
		return getImpl()->distributeDivOp(getOperation(), builder, resultType, value);
	}
};

namespace detail
{
	struct DerivativeInterfaceTraits
	{
		struct Concept
		{
			Concept() = default;
			Concept(const Concept& other) = default;
			Concept(Concept&& other) = default;
			Concept& operator=(Concept&& other) = default;
			virtual ~Concept() = default;
			Concept& operator=(const Concept& other) = default;

			virtual void derive(mlir::Operation* op, mlir::OpBuilder& builder, mlir::BlockAndValueMapping& derivatives) const = 0;
		};

		template <typename ConcreteOp>
		struct Model : public Concept
		{
			void derive(mlir::Operation* op, mlir::OpBuilder& builder, mlir::BlockAndValueMapping& derivatives) const final
			{
				return mlir::cast<ConcreteOp>(op).derive(builder, derivatives);
			}
		};

		template<typename ConcreteOp>
		class FallbackModel : public Concept
		{
			public:
			FallbackModel() = default;

			void derive(mlir::Operation* op, mlir::OpBuilder& builder, mlir::BlockAndValueMapping& derivatives) const final
			{
				return mlir::cast<ConcreteOp>(op).derive(builder, derivatives);
			}
		};
	};
}

class DerivativeInterface : public mlir::OpInterface<DerivativeInterface, detail::DerivativeInterfaceTraits>
{
	public:
	using OpInterface<DerivativeInterface, detail::DerivativeInterfaceTraits>::OpInterface;

	void derive(mlir::OpBuilder& builder, mlir::BlockAndValueMapping& derivatives)
	{
		return getImpl()->derive(getOperation(), builder, derivatives);
	}
};

namespace detail
{
	struct HeapAllocatorTraits
	{
		static llvm::StringRef getAutoFreeAttrName()
		{
			return "auto_free";
		}

		struct Concept
		{
			Concept() = default;
			Concept(const Concept& other) = default;
			Concept(Concept&& other) = default;
			Concept& operator=(Concept&& other) = default;
			virtual ~Concept() = default;
			Concept& operator=(const Concept& other) = default;

			static llvm::StringRef getAutoFreeAttrName()
			{
				return HeapAllocatorTraits::getAutoFreeAttrName();
			}

			virtual bool shouldBeFreed(mlir::Operation* op) const
			{
				llvm::StringRef attrName = getAutoFreeAttrName();
				return op->template getAttrOfType<mlir::BoolAttr>(attrName).getValue();
			}

			void setAsAutomaticallyFreed(mlir::Operation* op)
			{
				auto attr = mlir::BoolAttr::get(op->getContext(), true);
				op->setAttr(getAutoFreeAttrName(), attr);
			}

			void setAsManuallyFreed(mlir::Operation* op)
			{
				auto attr = mlir::BoolAttr::get(op->getContext(), false);
				op->setAttr(getAutoFreeAttrName(), attr);
			}
		};

		template <typename ConcreteOp>
		struct Model : public Concept
		{
			bool shouldBeFreed(mlir::Operation* op) const override
			{
				llvm::StringRef attrName = getAutoFreeAttrName();
				return op->template getAttrOfType<mlir::BoolAttr>(attrName).getValue();
			}
		};

		template<typename ConcreteOp>
		class FallbackModel : public Concept
		{
			public:
			FallbackModel() = default;

			bool shouldBeFreed(mlir::Operation* op) const override
			{
				return true;
			}
		};
	};
}

class HeapAllocator : public mlir::OpInterface<HeapAllocator, detail::HeapAllocatorTraits>
{
	public:
	using OpInterface<HeapAllocator, detail::HeapAllocatorTraits>::OpInterface;

	template <typename ConcreteOp>
	struct HeapAllocatorTrait : public mlir::OpInterface<HeapAllocator, detail::HeapAllocatorTraits>::Trait<ConcreteOp>
	{
		static llvm::StringRef getAutoFreeAttrName()
		{
			return detail::HeapAllocatorTraits::getAutoFreeAttrName();
		}

		bool shouldBeFreed()
		{
			mlir::Operation* op = (*static_cast<ConcreteOp*>(this)).getOperation();
			return mlir::cast<HeapAllocator>(op).shouldBeFreed();
		}

		void setAsAutomaticallyFreed()
		{
			mlir::Operation* op = (*static_cast<ConcreteOp*>(this)).getOperation();
			return mlir::cast<HeapAllocator>(op).setAsAutomaticallyFreed();
		}

		void setAsManuallyFreed()
		{
			mlir::Operation* op = (*static_cast<ConcreteOp*>(this)).getOperation();
			return mlir::cast<HeapAllocator>(op).setAsManuallyFreed();
		}
	};

	template <typename ConcreteOp>
	struct Trait : public HeapAllocatorTrait<ConcreteOp> {};

	static llvm::StringRef getAutoFreeAttrName()
	{
		return detail::HeapAllocatorTraits::getAutoFreeAttrName();
	}

	bool shouldBeFreed()
	{
		return getImpl()->shouldBeFreed(getOperation());
	}

	void setAsAutomaticallyFreed()
	{
		getImpl()->setAsAutomaticallyFreed(getOperation());
	}

	void setAsManuallyFreed()
	{
		getImpl()->setAsManuallyFreed(getOperation());
	}
};
