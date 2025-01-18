String configName = "debian-12"
String dockerfile = "debian-12.Dockerfile"
String checkName = "docker-dev-debug-image"

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
            def scmVars = checkout(scm)
            env.GIT_COMMIT = scmVars.GIT_COMMIT
            tag = sh(script: 'git describe --tags --abbrev=0', returnStdout: true).trim()
        }
    }

    String dockerRuntimeImageName = 'marco-compiler/marco-runtime-dev-debug-' + configName

    String dockerArgs =
        " --build-arg LLVM_PARALLEL_COMPILE_JOBS=${LLVM_PARALLEL_COMPILE_JOBS}" +
        " --build-arg LLVM_PARALLEL_LINK_JOBS=${LLVM_PARALLEL_LINK_JOBS}" +
        " --build-arg LLVM_BUILD_TYPE=Release" +
        " --build-arg LLVM_ENABLE_ASSERTIONS=ON" +
        " -f " + runtimeSrcPath + "/.jenkins/" + dockerfile +
        " " + runtimeSrcPath + "/.jenkins";

    publishChecks(name: checkName, status: 'IN_PROGRESS', summary: 'In progress')

    def dockerImage

    stage('Build') {
        dockerImage = docker.build(dockerRuntimeImageName + ":" + env.GIT_COMMIT[0..6], dockerArgs)
    }

    docker.withRegistry('https://ghcr.io', 'marco-ci') {
        stage('Publish') {
            dockerImage.push()

            if (tag != "") {
                dockerImage.push(tag)
            }

            dockerImage.push("latest")
        }
    }

    publishChecks(name: checkName, conclusion: 'SUCCESS', summary: 'Completed')
}
