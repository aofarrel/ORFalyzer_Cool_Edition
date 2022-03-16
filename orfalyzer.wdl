version 1.0

task orf {
	input {
		File inFile
		String outFile
	}

	command {
		touch ~{outFile}
		python3 /orfalyzer/orfalyzer.py ~{inFile} ~{outFile}
	}

	output {
		File orfs = "~{outFile}"
	}

	runtime {
        docker: "quay.io/aofarrel/orfalyzer:python3.9.10slim"
    }

    meta {
        author: "Ash O'Farrell"
        email: "aofarrel@ucsc.edu"
    }

}

workflow orfalyzer {
	input {
		File inFile
		String outFile = "ORFs.txt"
	}
	call orf {
		input: 
			inFile = inFile,
			outFile = outFile
	}
}