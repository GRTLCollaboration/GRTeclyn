# Visualising outputs

---

NOTE: THE FOLLOWING EXAMPLE HAS NOT YET BEEN PORTED... UNDER CONSTRUCTION...

---

AMReX outputs checkpoint and plot files in hdf5 format. There are several options for viewing and processing such files, but we generally use VisIt, ParaView or yt.

## Using Visit

Download Visit to your local machine from their [current releases](https://visit-dav.github.io/visit-website/releases-as-tables/).

For Mac and Windows there are installers, for Linux you should download the tar file, plus the "Visit Install Script" (in the bullets above the executable) and follow the instructions in "Visit Install Notes". The tar file for Ubuntu 14.04 seems to work on Ubuntu 16.04 too.

Assuming your hdf5 data is on a remote cluster, you have three options:

1. Download the hdf5 files to a local machine (or more likely onto an external hard drive connected to it, since the files are large) and run them directly there. (The command for copying files is scp).

2. Install (or module load) VisIt and run it in command line mode either directly on the login nodes or by submitting a batch job. Some example scripts and the appropriate run command for this can be found here:

https://github.com/GRTLCollaboration/Postprocessing_tools/tree/master/VisItTools

This is usually the best option for systems with a firewall preventing outgoing connections (e.g. Marenostrum, SupermucNG). One can generate png files or videos which can then be transferred to one's local machine using scp or by mounting the remote server using the sshfs command.

3. Run Visit remotely by downloading **the same version** (ie, 1.12.3 etc) of Visit on the cluster and setting up a remote profile in Visit on your local machine (see below). This has the advantage that you can keep the data on the cluster, and use its (probably more powerful) compute power, although some clusters don't like you to run Visit on the login nodes as it can clog up the system for other users, and may have dedicated nodes for visualisation. You should check this with your local cluster administrator. **(Note: The remote version should be the one with Mesa support for rendering without a display, otherwise you will have problems saving images and movies.)**

Assuming you choose option 3, the following information will help you get started:

### Setting up a remote host

To set up a remote host, launch VisIt on your local machine, then go to "Options->Host profiles". Click on "New Host" and configure it by setting:

* The Host nickname e.g. cosmos
* The remote hostname, e.g. cosmos.damtp.cam.ac.uk
* If you can run in parallel on the cluster nodes, set the max number of nodes and processors to use
* Path to visit installation (where you put it on the cluster), e.g. ~/visit
* Username (your username on the remote host), e.g. kclough
* Select "Tunnel data connections through ssh" and set the ssh command to "ssh -C"

VisIt can sometimes be "difficult" when it comes to getting it to remember the configuration. After you've configured the host, click "Apply", then click on the "i" button on the bottom-right of the main Visit window. Return to the host configuration dialog and hit "Export host" you should then see a confirmation that VisIt has saved it for you.

A useful tutorial on running VisIt in parallel in client server mode, where one needs to submit a batch job to reserve the compute nodes, can be found [here](https://www.youtube.com/watch?v=IPekLxRJxJ0). (Note that for this you need to download a version which includes parallel VisIt - ie redhat not ubuntu).

### Making a plot

VisIt has a GUI interface, so it is (sort of) intuitive. Opening a file should allow you to view the series of hdf5 outputs as a time series, without having to select any special options.

The most useful plots for our data are Pseudocolour plots, using the Operators->Slice operators to view a slice (adjust the intercept to the centre of the grid) and Operators->Elevate to make the plot 3D, but the things you can do with VisIt are pretty limitless.

There are VisIt tutorials which will help to discover all the functionality, see for example the ATPESC 2016 workshop which also has a youtube video:

http://www.visitusers.org/index.php?title=Short_Tutorial

See our tips for making good visualisations in [[Visualisation tips | Visualisation tips]]. Feel free to add to them!

There are a number of example scripts for processing GRChombo files using the VisIt command line here:

https://github.com/GRTLCollaboration/Postprocessing_tools/tree/master/VisItTools

## Using ParaView

[ParaView](https://www.paraview.org/) is an open-source and scalable
visualisation application similar to VisIt (it actually uses VisIt's Chombo
file reader underneath).

You can download a prebuilt version of ParaView for your local machine
(Windows/macOS/Linux) from [here](https://www.paraview.org/download/).
* For Windows and macOS, these are in the form of executables
  (exe and dmg respectively).
* For Linux, you can download and extract the tar file to a directory of your
  choosing. The application can be run by changing to the the `bin` subdirectory
  of the extracted folder and running the `paraview` executable, for example
  ```bash
  cd /path/to/ParaView-5.9.1-MPI-Linux-Python3.8-64bit/bin
  ./paraview
  ```
Before downloading the latest version, read below on the guide to remote
visualisation as you may need to download an older version to make it work.

### Remote visualisation

Assuming you have just performed a simulation on a remote HPC cluster and wish
to visualise the outputted HDF5 files, the best way to do this is to set up
ParaView for remote visualisation (client/server mode). In order to do this, you
will need to match the version on the cluster with your local version
so check what version is available via modules (assuming your system uses some
form of modules):
```
module avail paraview
```
and install the corresponding version locally. If ParaView is not installed,
ask the cluster administrator to install the latest version for you.

Since ParaView relies on open ports in order to be able to connect between the
client and server, and most HPC systems do not leave ports open for security
reasons, we will get around this by "tunnelling" the port with ssh.

To do this, follow the steps below
1. Open ParaView locally.
2. Click the 'Connect' icon (<img src="https://gitlab.kitware.com/paraview/paraview/-/raw/3762f9acc03c221b1755b137820855edb148edc5/Qt/Components/Resources/Icons/pqConnect.svg" width=20>)
   near the top left (or via the menu option File → Connect). Click 'Add Server'
   and set the fields to
    | Field | Value |
    | --- | --- |
    | Name | localhost |
    | Server Type | Client/Server |
    | Host | localhost |
    | Port | 11111 |

    Then click 'Configure' and set the 'Startup Type' field to 'Manual'.
    Click 'Save' and then click 'Close'
3. SSH into the remote cluster and load the relevant ParaView module. For versions
   5.10 or greater, the maximum
   number of threads that ParaView uses is controlled by the `VTK_SMP_MAX_THREADS`
   environment variable (similar to the `OMP_NUM_THREADS` environment variable for
   OpenMP programs such as GRChombo) so, in order to avoid accidentally using a
   very large number of threads, you should manually set this to a sensible value
   by running a command such as
   ```bash
   export VTK_SMP_MAX_THREADS=4
   ```
   Next, run the command
   ```
   pvserver
   ```
   which should print out something like
   ```
   Waiting for client...
   Connection URL: cs://<hostname>:<remote port>
   Accepting connection(s): <hostname>:<remote port>
   ```
   where `<remote port>` is usually something `11111`. Note that if loading a
   particularly large file, you might want to run this command in parallel
   with MPI in a job:
   ```
   mpiexec -n 8 pvserver
   ```
4. Set up the tunnel between the local `11111` port and `<remote port>` with
   the command
   ```bash
   ssh -NL 11111:localhost:<remote port> <username>@<hostname>
   ```
   It may ask you to authenticate in the usual way and then look as though it
   has 'hung' (i.e. no prompt). If you are running `pvserver` in a job, and
   you cannot SSH directly to the compute node on which the job is running
   (as is likely the case),
   you will need to tunnel the port via the login node which can be done with
   a command such as
   ```bash
   ssh -L 11111:localhost:11111 <username>@<login node hostname> ssh -4 -L 11111:localhost:<remote port> -N <compute node hostname>
   ```
5. Click the 'Connect' icon again and choose the server we configured in step 2
   called `localhost`. It should then connect to the remote cluster and the
   output from your SSH session in step 3 will have the extra line
   ```
   Client connected.
   ```
   Using the menu options 'File → Open', you should be able to browse the remote
   filesystem and select your HDF5 files.

Note that you only need to do step 2 once. To run remote visualisation another
time, simply repeat steps 1 and 3-5.

### Remote visualisation (reverse connection)
If you are having problems with the above method for remote visualisation with
ParaView, there is an alternative in the form of *reverse connections* where
the remote server connects to the local client (rather than the above where the
local client connects to the server). If you are trying to use Catalyst Live
with the ParaView Catalyst insitu instrumentation, this also uses a reverse
connection so you will need to follow similar steps.

As for the conventional client/server
mode, you will need to have a version of ParaView installed on the remote system
and the same version installed locally.

To set up remote visualisation with reverse connections, follow the steps below
1. Open ParaView locally
2. Click the 'Connect' icon (<img src="https://gitlab.kitware.com/paraview/paraview/-/raw/3762f9acc03c221b1755b137820855edb148edc5/Qt/Components/Resources/Icons/pqConnect.svg" width=20>)
   near the top left (or via the menu option File → Connect). Click 'Add Server'
   and set the fields to
    | Field | Value |
    | --- | --- |
    | Name | localhost (reverse connection) |
    | Server Type | Client/Server (reverse connection) |
    | Port | 11111 |

   Then click 'Configure' and set the 'Startup Type' field to 'Manual'.
   Click 'Save'.
3. Connect to the server we have just configured by selecting it and then
   clicking 'Connect'. A dialog box will appear which says
   > Establishing connection to 'localhost (reverse connection)'
   > Waiting for server to connect.
4. Set up the tunnel between the local `11111` port and the remote port `xxxxx`
   (this can usually be set to `11111` but we will leave it generic in the
   following instructions) with the command
   ```bash
   ssh -NR 11111:localhost:xxxxx username@hostname
   ```
   It may ask for authentication and then look as though it has 'hung' (i.e. no
   prompt) If you wish to run the ParaView server in a job, and you cannot SSH
   directly to the compute node on which the job is running, you will need to
   tunnel the port via the login node which can be done with a command such as
   ```bash
   ssh -R 11111:localhost:yyyyy <username>@<login node hostname> ssh -R yyyyy:localhost:xxxxx -N <compute node hostname>
   ```
   Note that these commands are virtually identical to the ones for the
   conventional client/server tunnelling but the `-L` flag has changed to `-R`.
5. SSH into the remote cluster and load the relevant ParaView module. For
   versions 5.10 or greater, the maximum number of threads that ParaView uses is
   controlled by the `VTK_SMP_MAX_THREADS` environment variable (similar to the
   `OMP_NUM_THREADS` environment variable for OpenMP programs such as GRChombo)
   so, in order to avoid accidentally using a very large number of threads, you
   should manually set this to a sensible value by running a command such as
   ```bash
   export VTK_SMP_MAX_THREADS=4
   ```
   Next, run the command
   ```
   pvserver --reverse-connection --server-port=xxxxx
   ```
   or alternatively in parallel (preferably in a jobscript)
   ```
   mpiexec -n 4 pvserver --reverse-connection --server-port=xxxxx
   ```
   Assuming everything has worked, you should get the following output
   ```
   Connecting to client (reverse connection requested)...
   Connection URL: csrc://localhost:xxxxx
   Client connected.
   ```
   Using the menu options 'File → Open', you should be able to browse the remote
   filesystem and select your HDF5 files.


### Documentation and Tutorials

The ParaView user and reference guide can be found
[here](https://docs.paraview.org/en/latest/). Make sure to select the correct
version in the bottom left.

There are also some tutorials that are linked to from the main ParaView website
[here](https://www.paraview.org/tutorials/).

As for VisIt above, there is also a YouTube video for a presentation given at
ATPESC [here](https://youtu.be/sXY72e3Ce4g).


## Using yt

yt is an alternative python based visualisation software which is very good for processing and analysing data. Instructions on downloading and using yt can be found here:

http://yt-project.org/doc/index.html

There are a number of example scripts for processing GRChombo files here:

https://github.com/GRChombo/Postprocessing_tools/tree/master/YTAnalysisTools

### Basics

The following script shows an example of the most basic commands:

```
import yt

filename = "/your/path/BBH_000100.3d.hdf5"       # path to the checkpoint/plot file.
ds = yt.load(filename)                           # yt.load() automatically detects that is a Chombo file and loads it.

L, _, _ = ds.domain_with                         # extract the size of the grid

# Loading data
data_flat = ds.r[:,:,:]                          # creates a dict with flat data (1D-array),
                                                 #it ignores duplicated datapoints from coarser levels.
data_grid = ds.r[::120j,::120j,::120j]           # creates a dict with grid data (3d-array with 120 points per side),
                                                 ## using 0th-interpolation order ('nearest') from the finest level.
data_grid = ds.r[::120j,::120j, L/2]             # creates a dict with grid data (2d-array with 120 points per side),
                                                 ## using 0th-interpolation order ('nearest') from the finest level.


print('shape flat: ',  data_flat["K"].shape)     # Here it has been used the var "K" as an example
print('shape grid: ',  data_grid["K"].shape)
print('shape slice: ',  data_slice["K"].shape)
# Output:
# shape flat:  (15489664,)
# shape grid:  (120, 120, 120)
# shape slice:  (120, 120)

```

The command `ds.r[:,:,:]` creates a python-dictionary that contains the outputted variables  (e.g. "K", "chi", etc) and other useful grid-variables  (e.g. "x", "y", ..., "dx", ..., etc).
The following script shows an example of how to compute the averaged quantities of a variable of interest.

```
import numpy as np

# Loading data
ds = yt.load("/path/your/file/???.hdf5")
dd = ds.r[:,:,:]                                         # Load the dict containing the flat array data

gridcell_volume = dd['dx']**3
physical_cell_volume =  dd['dx']**3*dd['chi']**(-1.5)    # physical cell volume taking into account the conformal factor "chi".
total_volume = np.sum(physical_cell_volume)

average_K = np.sum(dd['K'] * physical_cell_volume)/total_volume

```

yt contain many additional functionalities for both analysis and plotting. Have a look at their documentation for an in-depth description (link above).