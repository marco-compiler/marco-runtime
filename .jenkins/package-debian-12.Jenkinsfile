String configName = "debian-12"
String dockerfile = "debian-12.Dockerfile"

node {
    agent {
        label 'x86_64-linux'
    }

    String localWorkspace = "${WORKSPACE}/" + configName

    String srcPath = localWorkspace + "/src"
    String buildPath = localWorkspace + "/build"
    String installPath = localWorkspace + "/install"

    String runtimeSrcPath = srcPath + "/marco-runtime"
    String runtimeBuildPath = buildPath + "/marco-runtime"
    String runtimeInstallPath = installPath + "/marco-runtime"

    stage("Checkout") {
        dir(runtimeSrcPath) {
            checkout(scm)
        }
    }

    String dockerRuntimeImageName = 'marco-compiler/marco-runtime-package-' + configName

    String dockerArgs =
        " --build-arg LLVM_PARALLEL_COMPILE_JOBS=${LLVM_PARALLEL_COMPILE_JOBS}" +
        " --build-arg LLVM_PARALLEL_LINK_JOBS=${LLVM_PARALLEL_LINK_JOBS}" +
        " --build-arg LLVM_BUILD_TYPE=Release" +
        " --build-arg LLVM_ENABLE_ASSERTIONS=OFF" +
        " -f " + runtimeSrcPath + "/.jenkins/" + dockerfile +
        " " + runtimeSrcPath + "/.jenkins";

    stage('Docker image') {
        dockerImage = docker.build(dockerRuntimeImageName + ':latest', dockerArgs)
    }

    dockerImage.inside() {
        stage("OS information") {
            sh "cat /etc/os-release"
        }

        stage('Configure') {
            cmake arguments: "-S " + runtimeSrcPath + " -B " + runtimeBuildPath + " -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_LINKER_TYPE=MOLD -DCMAKE_INSTALL_PREFIX=" + runtimeInstallPath, installation: 'InSearchPath', label: 'Configure'
        }

        stage('Install') {
            cmake arguments: "--build " + runtimeBuildPath + " --target install", installation: 'InSearchPath', label: 'Install'
        }

        stage('Package') {
            sh "chmod +x " + runtimeSrcPath + "/.jenkins/package/" + configName + "/build.sh"
            sh runtimeSrcPath + "/.jenkins/package/" + configName + "/build.sh " + runtimeSrcPath + " " + runtimeInstallPath

            sshPublisher(publishers: [sshPublisherDesc(configName: 'marco-package', transfers: [sshTransfer(cleanRemote: false, excludes: '', execCommand: '', execTimeout: 120000, flatten: false, makeEmptyDirs: false, noDefaultExcludes: false, patternSeparator: '[, ]+', remoteDirectory: configName + "/amd64", remoteDirectorySDF: false, sourceFiles: '*.deb')], usePromotionTimestamp: false, useWorkspaceInPromotion: false, verbose: false)])
        }
    }
}
