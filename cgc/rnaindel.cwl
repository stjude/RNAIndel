{
    "class": "CommandLineTool",
    "cwlVersion": "v1.2",
    "$namespaces": {
        "sbg": "https://sevenbridges.com"
    },
    "baseCommand": [],
    "inputs": [
        {
            "id": "subcommand",
            "type": {
                "type": "enum",
                "symbols": [
                    "PredictIndels",
                    "CalculateFeatures",
                    "Train",
                    "CountOccurrence"
                ],
                "name": "subcommand"
            },
            "inputBinding": {
                "shellQuote": false,
                "position": 0
            },
            "default": "PredictIndels"
        },
        {
            "id": "input",
            "type": "File",
            "inputBinding": {
                "prefix": "-i",
                "shellQuote": false,
                "position": 1
            },
            "doc": "STAR-mapped BAM file to analyze",
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
            "loadListing": "deep_listing",
            "id": "data_dir",
            "type": "Directory",
            "inputBinding": {
                "prefix": "-d",
                "shellQuote": false,
                "position": 2
            },
            "doc": "Data directory containing trained models and databases. Can be obtained from http://ftp.stjude.org/pub/software/RNAIndel/data_dir_grch38.v3.tar.gz (GRCh38) or http://ftp.stjude.org/pub/software/RNAIndel/data_dir_grch37.v3.tar.gz (GRCh37)"
        },
        {
            "id": "refdata",
            "type": "File",
            "inputBinding": {
                "prefix": "-r",
                "shellQuote": false,
                "position": 3
            },
            "doc": "Reference genome in FASTA format",
            "sbg:fileTypes": "FA",
            "secondaryFiles": [
                {
                    "pattern": ".fai",
                    "required": true
                }
            ]
        },
        {
            "id": "p",
            "type": "int?",
            "inputBinding": {
                "prefix": "-p",
                "shellQuote": false,
                "position": 5
            },
            "label": "number of cores",
            "default": 8
        }
    ],
    "outputs": [
        {
            "id": "predicted_indels",
            "type": "File",
            "outputBinding": {
                "glob": "*.vcf"
            }
        }
    ],
    "doc": "RNAIndel calls coding indels from tumor RNA-Seq data and classifies them as somatic, germline, and artifactual. RNAIndel supports GRCh38 and 37.\n\n## Inputs\n* **BAM** - STAR-mapped BAM file\n* **Fasta** - Reference genome in FASTA format\n* **Reference** - Trained data models and databases. Can be obtained from http://ftp.stjude.org/pub/software/RNAIndel/data_dir_grch38.v3.tar.gz (GRCh38) or http://ftp.stjude.org/pub/software/RNAIndel/data_dir_grch37.v3.tar.gz (GRCh37)\n\n## Outputs\n* **Indel callset** - RNAIndel called indels",
    "label": "rnaindel2",
    "arguments": [
        {
            "prefix": "-o",
            "shellQuote": false,
            "position": 5,
            "valueFrom": "$(inputs.input.nameroot).vcf"
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
            "class": "ResourceRequirement",
            "ramMin": 80000,
            "coresMin": 8
        },
        {
            "class": "DockerRequirement",
            "dockerPull": "cgc-images.sbgenomics.com/stjude/rnaindel:latest"
        },
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "hints": [
        {
            "class": "sbg:AWSInstanceType",
            "value": "r4.4xlarge;ebs-gp2;1024"
        }
    ],
    "sbg:links": [
        {
            "id": "https://github.com/stjude/RNAIndel",
            "label": "Source Code"
        },
        {
            "id": "https://doi.org/10.1093/bioinformatics/btz753",
            "label": "Publication"
        }
    ],
    "sbg:appVersion": [
        "v1.2"
    ],
    "sbg:wrapperLicense": "Apache 2.0 License",
    "sbg:license": "Apache 2.0 License",
    "sbg:categories": [
        "RNA-Seq",
        "Variant Calling"
    ]
}
