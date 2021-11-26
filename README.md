# SPAtial DEconvolution Benchmark
This is a repository for benchmarking spatial deconvolution and mapping tools.

You can either run the pipeline with the gold/silver/bronze standards or run it with your own data. In either case, you need to [install NextFlow](https://www.nextflow.io/docs/latest/getstarted.html), and either [Docker](https://docs.docker.com/get-docker/) or [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html).

In the future, you will be able to run this pipeline directly from NextFlow. From now you must clone the repository, `cd spade-benchmark` and run scripts from there.

### Reproducing the pipeline with standards



### Platforms
The workflow has been tested on two platforms:
- NextFlow 21.04.3 on the Windows Subsystem for Linux (WSL2, Ubuntu 20.04), with Docker Desktop for Windows 4.1.0
- NextFlow 20.10.3 on CentOS Linux, with Singularity 3.8.1
