String configName = "debian-12"
String dockerfile = "prod-debian-12.Dockerfile"
String checkName = "docker-prod-image"

publishChecks(name: checkName, status: 'QUEUED', summary: 'Queued')

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

    String dockerRuntimeImageName = 'marco-compiler/marco-runtime-prod-' + configName

    String dockerArgs =
        " --build-arg LLVM_PARALLEL_COMPILE_JOBS=${LLVM_PARALLEL_COMPILE_JOBS}" +
        " --build-arg LLVM_PARALLEL_LINK_JOBS=${LLVM_PARALLEL_LINK_JOBS}" +
        " -f " + runtimeSrcPath + "/.jenkins/" + dockerfile +
        " " + runtimeSrcPath + "/.jenkins";

    publishChecks(name: checkName, status: 'IN_PROGRESS', summary: 'In progress')

    def dockerImage

    stage('Build') {
        dockerImage = docker.build(dockerRuntimeImageName + ':latest', dockerArgs)
    }

    docker.withRegistry('https://ghcr.io', 'marco-ci') {
        stage('Publish') {
            dockerImage.push()
        }
    }

    publishChecks(name: checkName, conclusion: 'SUCCESS', summary: 'Completed')
}
