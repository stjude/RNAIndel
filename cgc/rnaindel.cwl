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
                    "SetUp",
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
        },
        {
            "id": "tumor",
            "type": "File?",
            "inputBinding": {
                "prefix": "-t",
                "shellQuote": false,
                "position": 6
            },
            "label": "Tumor DNA-Seq BAM file for cross-platform check (default: None)"
        },
        {
            "id": "normal",
            "type": "File?",
            "inputBinding": {
                "prefix": "-n",
                "shellQuote": false,
                "position": 7
            },
            "label": "Normal DNA-Seq BAM file for cross-platform check (default: None)"
        },
        {
            "id": "vcf",
            "type": "File?",
            "inputBinding": {
                "prefix": "-v",
                "shellQuote": false,
                "position": 8
            },
            "label": "VCF file from external caller. Supply as vcf.gz + index",
            "secondaryFiles": [
                {
                    "pattern": ".tbi",
                    "required": true
                }
            ]
        },
        {
            "id": "mapq",
            "type": "int?",
            "inputBinding": {
                "prefix": "-q",
                "shellQuote": false,
                "position": 9
            },
            "label": "STAR mapping quality MAPQ for unique mappers (default: 255)",
            "default": 255
        },
        {
            "id": "heap",
            "type": "string?",
            "inputBinding": {
                "prefix": "-m",
                "shellQuote": false,
                "position": 10
            },
            "label": "maximum heap space (default: 6000m)",
            "default": "6000m"
        },
        {
            "id": "pon",
            "type": "File?",
            "inputBinding": {
                "prefix": "--pon",
                "shellQuote": false,
                "position": 11
            },
            "label": "User defined panel of normals to refine somatic predictions. Supply as vcf.gz + index",
            "secondaryFiles": [
                {
                    "pattern": ".tbi",
                    "required": true
                }
            ]
        },
        {
            "id": "region",
            "type": "string?",
            "inputBinding": {
                "prefix": "--region",
                "shellQuote": false,
                "position": 12
            },
            "label": "specify region for target analysis: chrN:start-stop (default: None)"
        },
        {
            "id": "sensitive_mode",
            "type": "boolean?",
            "default": False,
            "inputBinding": {
                "prefix": "--deactivate-sensitive-mode",
                "shellQuote": false,
                "position": 13
            },
            "label": "Disable additional realignments for soft-clipped reads (default: False)"
        },
        {
            "id": "safety_mode",
            "type": "boolean?",
            "default": False,
            "inputBinding": {
                "prefix": "--safety-mode",
                "shellQuote": false,
                "position": 14
            },
            "label": "Deactivate parallelism at realignment step. may be required to run with -p > 1 on some platforms. (default: False)"
        },
        {
            "id": "skip_homopolyer_outlier_analysis",
            "type": "boolean?",
            "default": False,
            "inputBinding": {
                "prefix": "--skip-homopolyer-outlier-analysis",
                "shellQuote": false,
                "position": 15
            },
            "label": "No outlier analysis for homopolymer indels (repeat > 4) performed if set. (default: False)"
        },
        {
           "id": "include_all_external_calls",
            "type": "boolean?",
            "default": False,
            "inputBinding": {
                "prefix": "--include-all-external-calls",
                "shellQuote": false,
                "position": 16
            },
            "label": "Set to include all indels in VCF file supplied by -v. (default: False. Use only calls with PASS in FILTER)"
        }
    ],
    "outputs": [
        {
            "id": "predicted_indels",
            "type": "File",
            "outputBinding": {
                "glob": "*.vcf.gz"
            },
            "secondaryFiles": [
                {
                    "pattern": ".tbi",
                    "required": true
                }
            ]
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
