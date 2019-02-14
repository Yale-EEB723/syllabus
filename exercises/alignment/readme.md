
# Running the exercise

The script in this directory (`alignment.sh`) should run on its own as long as you create a directory for it called `/data/eeb723-seqaln`. However, we recommend that you open it in a separate window then copy and paste commands into an interactive shell one by one. Below are instructions for creating a shell with Docker on your computer or Singularity on the Farnam cluster.

## Docker

### Windows

Open a new command prompt by typing <kbd><img src=http://i.stack.imgur.com/B8Zit.png></kbd>+<kbd>R</kbd>, then enter `cmd` into the prompt and click OK


We need a folder to keep files generated from our analyses on your computer. Enter the following into your command prompt:

``` cmd
md %USERPROFILE%\Desktop\eeb723-seqaln
```

Next, start up the class environment in Docker:

```
docker run --rm -ti -v %USERPROFILE%\Desktop\eeb723-seqaln:/data/eeb723-seqaln eeb723/course_docker bash
```


### macOS/linux

Open a new Terminal window (on a mac start up Terminal.app), then enter the following to set up a directory to keep files generated from our analyses on your computer.

``` bash
mkdir -p ~/Desktop/eeb723-seqaln
docker run --rm -ti -v ~/Desktop/eeb723-seqaln:/data/eeb723-seqaln eeb723/course_docker bash
```


## Singularity on Farnam

Connect to Farnam (on windows use [mobaxterm](http://docs.ycrc.yale.edu/clusters-at-yale/access/#connect-from-windows))

``` bash
ssh netid@farnam.hpc.yale.edu
```

From the login node, allocate an interactive job with `srun`.

``` bash
# interactive session with 2 cores and default 10GiB of RAM
srun -A eeb723 --pty -p interactive -c 2 bash
```

Choose a working directory and band-mount it to your container

``` bash
# if you prefer a different directory use that instead
mkdir -p /gpfs/ysm/project/eeb723/${USER}/eeb723-seqaln
singularity shell --shell /bin/bash -B /gpfs/ysm/project/eeb723/${USER}/eeb723-seqaln:/data/eeb723-seqaln docker://eeb723/course_docker
```
 
