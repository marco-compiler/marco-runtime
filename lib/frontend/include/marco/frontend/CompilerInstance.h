#ifndef MARCO_FRONTEND_COMPILERINSTANCE_H
#define MARCO_FRONTEND_COMPILERINSTANCE_H

//#include "flang/Parser/parsing.h"
//#include "flang/Parser/provenance.h"
//#include "flang/Semantics/semantics.h"

#include <clang/Basic/DiagnosticOptions.h>
#include <list>
#include <llvm/Support/raw_ostream.h>

#include "CompilerInvocation.h"
#include "CompilerAction.h"
//#include "FrontendAction.h"

namespace marco::frontend
{
  /// Helper class for managing a single instance of the Flang compiler.
  ///
  /// This class serves two purposes:
  ///  (1) It manages the various objects which are necessary to run the compiler
  ///  (2) It provides utility routines for constructing and manipulating the
  ///      common Flang objects.
  ///
  /// The compiler instance generally owns the instance of all the objects that it
  /// manages. However, clients can still share objects by manually setting the
  /// object and retaking ownership prior to destroying the CompilerInstance.
  ///
  /// The compiler instance is intended to simplify clients, but not to lock them
  /// in to the compiler instance for everything. When possible, utility functions
  /// come in two forms; a short form that reuses the CompilerInstance objects,
  /// and a long form that takes explicit instances of any required objects.
  class CompilerInstance
  {
      /// The options used in this compiler instance.
      std::shared_ptr<CompilerInvocation> invocation_;

      /// Flang file manager.
      /*
      std::shared_ptr<Fortran::parser::AllSources> allSources_;

      std::shared_ptr<Fortran::parser::AllCookedSources> allCookedSources_;

      std::shared_ptr<Fortran::parser::Parsing> parsing_;

      std::unique_ptr<Fortran::semantics::Semantics> semantics_;
       */

      /// The stream for diagnostics from Semantics
      llvm::raw_ostream* semaOutputStream_ = &llvm::errs();

      /// The stream for diagnostics from Semantics if owned, otherwise nullptr.
      std::unique_ptr<llvm::raw_ostream> ownedSemaOutputStream_;

      /// The diagnostics engine instance.
      llvm::IntrusiveRefCntPtr<clang::DiagnosticsEngine> diagnostics_;

      /// Holds information about the output file.
      struct OutputFile
      {
        std::string filename_;

        OutputFile(std::string inputFilename)
            : filename_(std::move(inputFilename)) {}
      };

      /// The list of active output files.
      std::list<OutputFile> outputFiles_;

      /// Holds the output stream provided by the user. Normally, users of
      /// CompilerInstance will call CreateOutputFile to obtain/create an output
      /// stream. If they want to provide their own output stream, this field will
      /// facilitate this. It is optional and will normally be just a nullptr.
      std::unique_ptr<llvm::raw_pwrite_stream> outputStream_;

    public:
      explicit CompilerInstance();

      ~CompilerInstance();

      CompilerInvocation& invocation()
      {
        assert(invocation_ && "Compiler instance has no invocation!");
        return *invocation_;
      };

      /// Replace the current invocation.
      void set_invocation(std::shared_ptr<CompilerInvocation> value);

      /// Return the current allSources.
      /*
      Fortran::parser::AllSources& allSources() const { return *allSources_; }

      bool HasAllSources() const { return allSources_ != nullptr; }

      parser::AllCookedSources& allCookedSources()
      {
        assert(allCookedSources_ && "Compiler instance has no AllCookedSources!");
        return *allCookedSources_;
      };

      /// }
      /// @name Parser Operations
      /// {

      /// Return parsing to be used by Actions.
      Fortran::parser::Parsing& parsing() const { return *parsing_; }
       */

      /// Replace the current stream for verbose output.
      void set_semaOutputStream(llvm::raw_ostream& Value);

      /// Replace the current stream for verbose output.
      void set_semaOutputStream(std::unique_ptr<llvm::raw_ostream> Value);

      /// Get the current stream for verbose output.
      llvm::raw_ostream& semaOutputStream() { return *semaOutputStream_; }

      bool ExecuteAction(CompilerAction& act);

      clang::DiagnosticOptions& GetDiagnosticOpts()
      {
        return invocation_->GetDiagnosticOpts();
      }

      const clang::DiagnosticOptions& GetDiagnosticOpts() const
      {
        return invocation_->GetDiagnosticOpts();
      }

      FrontendOptions& frontendOpts() { return invocation_->frontendOpts(); }

      const FrontendOptions& frontendOpts() const
      {
        return invocation_->frontendOpts();
      }

      /// Get the current diagnostics engine.
      clang::DiagnosticsEngine& diagnostics() const
      {
        assert(diagnostics_ && "Compiler instance has no diagnostics!");
        return *diagnostics_;
      }

      clang::DiagnosticConsumer& GetDiagnosticClient() const
      {
        assert(diagnostics_ && diagnostics_->getClient() &&
            "Compiler instance has no diagnostic client!");
        return *diagnostics_->getClient();
      }

      /// Clear the output file list.
      void ClearOutputFiles(bool eraseFiles);

      /// Create the default output file (based on the invocation's options) and
      /// add it to the list of tracked output files. If the name of the output
      /// file is not provided, it will be derived from the input file.
      ///
      /// \param binary     The mode to open the file in.
      /// \param baseInput  If the invocation contains no output file name (i.e.
      ///                   outputFile in FrontendOptions is empty), the input path
      ///                   name to use for deriving the output path.
      /// \param extension  The extension to use for output names derived from
      ///                   \p baseInput.
      /// \return           Null on error, ostream for the output file otherwise
      std::unique_ptr<llvm::raw_pwrite_stream> CreateDefaultOutputFile(
          bool binary = true, llvm::StringRef baseInput = "",
          llvm::StringRef extension = "");

    private:
      /// Create a new output file
      ///
      /// \param outputPath   The path to the output file.
      /// \param binary       The mode to open the file in.
      /// \return             Null on error, ostream for the output file otherwise
      llvm::Expected<std::unique_ptr<llvm::raw_pwrite_stream>> CreateOutputFileImpl(
          llvm::StringRef outputPath, bool binary);

    public:
      /// }
      /// @name Construction Utility Methods
      /// {

      /// Create a DiagnosticsEngine object
      ///
      /// If no diagnostic client is provided, this method creates a
      /// DiagnosticConsumer that is owned by the returned diagnostic object. If
      /// using directly the caller is responsible for releasing the returned
      /// DiagnosticsEngine's client eventually.
      ///
      /// \param opts - The diagnostic options; note that the created text
      /// diagnostic object contains a reference to these options.
      ///
      /// \param client - If non-NULL, a diagnostic client that will be attached to
      /// (and optionally, depending on /p shouldOwnClient, owned by) the returned
      /// DiagnosticsEngine object.
      ///
      /// \return The new object on success, or null on failure.
      static clang::IntrusiveRefCntPtr<clang::DiagnosticsEngine> CreateDiagnostics(
          clang::DiagnosticOptions* opts,
          clang::DiagnosticConsumer* client = nullptr, bool shouldOwnClient = true);

      void CreateDiagnostics(
          clang::DiagnosticConsumer* client = nullptr, bool shouldOwnClient = true);

      /// }
      /// @name Output Stream Methods
      /// {
      void set_outputStream(std::unique_ptr<llvm::raw_pwrite_stream> outStream)
      {
        outputStream_ = std::move(outStream);
      }

      bool IsOutputStreamNull() { return (outputStream_ == nullptr); }

      // Allow the frontend compiler to write in the output stream.
      void WriteOutputStream(const std::string& message)
      {
        *outputStream_ << message;
      }
  };
}

#endif // MARCO_FRONTEND_COMPILERINSTANCE_H
