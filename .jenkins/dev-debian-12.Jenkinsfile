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
            def scmVars = checkout(scm)
            env.GIT_COMMIT = scmVars.GIT_COMMIT
        }
    }

    String dockerRuntimeImageName = 'marco-compiler/marco-runtime-dev-' + configName

    String dockerArgs =
        " --build-arg LLVM_PARALLEL_COMPILE_JOBS=${LLVM_PARALLEL_COMPILE_JOBS}" +
        " --build-arg LLVM_PARALLEL_LINK_JOBS=${LLVM_PARALLEL_LINK_JOBS}" +
        " --build-arg LLVM_BUILD_TYPE=Release" +
        " --build-arg LLVM_ENABLE_ASSERTIONS=ON" +
        " -f " + runtimeSrcPath + "/.jenkins/" + dockerfile +
        " " + runtimeSrcPath + "/.jenkins";

    def dockerImage

    stage("Docker image") {
        dockerImage = docker.build(dockerRuntimeImageName + ':' + env.GIT_COMMIT, dockerArgs)
    }

    dockerImage.inside() {
        stage("OS information") {
            sh "cat /etc/os-release"
        }

        stage('Configure') {
            cmake arguments: "-S " + runtimeSrcPath + " -B " + runtimeBuildPath + " -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_LINKER_TYPE=MOLD -DCMAKE_INSTALL_PREFIX=" + runtimeInstallPath, installation: 'InSearchPath', label: 'Configure'
        }

        stage('Build') {
            cmake arguments: "--build " + runtimeBuildPath, installation: 'InSearchPath', label: 'Build'
        }

        stage('Unit test') {
            cmake arguments: "--build " + runtimeBuildPath + " --target test", installation: 'InSearchPath', label: 'Unit tests'
        }

        stage('Install') {
            cmake arguments: "--build " + runtimeBuildPath + " --target install", installation: 'InSearchPath', label: 'Install'
        }
    }
}
