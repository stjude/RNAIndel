{
    "class": "CommandLineTool",
    "cwlVersion": "v1.1",
    "id": "rnaindel",
    "baseCommand": [],
    "inputs": [
        {
            "id": "subcommand",
            "type": {
                "type": "enum",
                "symbols": [
                    "analysis",
                    "feature",
                    "nonsomatic",
                    "reclassification",
                    "recurrence",
                    "training"
                ],
                "name": "subcommand"
            },
            "inputBinding": {
                "shellQuote": false,
                "position": 0
            }
        },
        {
            "id": "input_bam",
            "type": "File",
            "inputBinding": {
                "prefix": "-i",
                "shellQuote": false,
                "position": 1
            },
            "label": "Input BAM file",
            "sbg:fileTypes": "BAM",
            "secondaryFiles": [
                {
                    "pattern": ".bai?"
                },
                {
                    "pattern": "$(self.nameroot).bai?"
                }
            ]
        },
        {
            "id": "fasta",
            "type": "File",
            "inputBinding": {
                "prefix": "-r",
                "shellQuote": false,
                "position": 3
            },
            "label": "Input FASTA",
            "sbg:fileTypes": "FA",
            "secondaryFiles": [
                {
                    "pattern": ".fai"
                }
            ]
        },
        {
            "loadListing": "deep_listing",
            "id": "data",
            "type": "Directory",
            "inputBinding": {
                "prefix": "-d",
                "shellQuote": false,
                "position": 4
            },
            "label": "Data Directory"
        },
        {
            "id": "user_caller_vcf",
            "type": "File?",
            "inputBinding": {
                "prefix": "-v",
                "shellQuote": false,
                "position": 10
            },
            "sbg:fileTypes": "VCF"
        },
        {
            "id": "star_mapping",
            "type": "int?",
            "inputBinding": {
                "prefix": "-q",
                "shellQuote": false,
                "position": 10
            }
        },
        {
            "id": "num_cores",
            "type": "int?",
            "inputBinding": {
                "prefix": "-p",
                "shellQuote": false,
                "position": 10
            }
        },
        {
            "id": "max_heap_space",
            "type": "string?",
            "inputBinding": {
                "prefix": "-m",
                "shellQuote": false,
                "position": 10
            }
        },
        {
            "id": "log_directory",
            "type": "string?",
            "inputBinding": {
                "prefix": "-l",
                "shellQuote": false,
                "position": 10
            }
        },
        {
            "id": "nonsomatic_indels_vcf",
            "type": "File?",
            "inputBinding": {
                "prefix": "-n",
                "shellQuote": false,
                "position": 10
            },
            "sbg:fileTypes": "VCF"
        },
        {
            "id": "germline_indels_vcf",
            "type": "File?",
            "inputBinding": {
                "prefix": "-g",
                "shellQuote": false,
                "position": 10
            },
            "sbg:fileTypes": "VCF"
        },
        {
            "id": "input",
            "type": "string?",
            "inputBinding": {
                "prefix": "--region",
                "shellQuote": false,
                "position": 10
            }
        }
    ],
    "outputs": [
        {
            "id": "#output_file",
            "type": "File",
            "outputBinding": {
                "glob": "*.vcf"
            }
        }
    ],
    "label": "rnaindel",
    "arguments": [
        {
            "prefix": "-o",
            "shellQuote": false,
            "position": 101,
            "valueFrom": "$(inputs.input_bam.nameroot).vcf"
        }
    ],
    "requirements": [
        {
            "class": "ShellCommandRequirement"
        },
        {
            "class": "LoadListingRequirement"
        },
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "hints": [
        {
            "class": "DockerRequirement",
            "dockerPull": "ghcr.io/stjude/rnaindel:latest"
        }
    ]
}
