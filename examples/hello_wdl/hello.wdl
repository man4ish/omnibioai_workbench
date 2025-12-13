version 1.0

workflow hello_wdl {
    call hello
}

task hello {
    command {
        echo "Hello from OmnibioAI WDL workflow"
    }
    output {
        String message = read_string(stdout())
    }
}

