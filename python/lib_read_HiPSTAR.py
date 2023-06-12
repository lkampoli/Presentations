"""
Library created by Javi Otero on 10/07/2017.

This library contains conveninent functions to read-in multi-block data from HiPSTAR
into python. This includes FLOW files, Subspaces, Sfiles and STAT_cont.bin files.
This is achieved by storing the blocks as a list of objects, allowing for different
dimensions across blocks. The object structure can be seen in class new_block, where
the data (either the grid coordinates or flowfield) is stored in self.data. As well,
this library allows to read only a subset of the data contained into one file,
which comes in very handy when dealing with huge datasets.

  ** Classes:

    - new_block: Class to store single-block data.


  ** Reading functions: They return the list of objects (instances of the class
                        new_block) where the multi-block data is stored.

    - To read FLOW files and Subspaces use:
       --> get_grid
       --> get_data

    - To read STAT_cont.bin file, the STAT.xmf wrapper needs to be created first.
      Then use:
       --> get_grid_stat_xmf
       --> get_data_stat_xmf

    - To read Sfiles use:
       --> read_Sfile
       --> averaged_Sfile
       --> averaged_Sfile_convergence

  ** Data processing functions: They carry out operations in the datasets created
                                by the functions described above.

    - To process Sfiles use:
       --> process_Sfile

    - To compute 3D mult-iblock metrics use:
       --> grid_metrics

  ** Writing functions:

    - To readAndWrite Subspaces(to be used for extracting a slice or subvolume of dataset):
       --> getAndWrite_grid
       --> getAndWrite_data

    - To write Sfiles use:
       --> write_Sfile
       --> merge_Sfile

    - To write STAT_cont files use:
       --> write_STAT_cont


Additionally, this library also contains functions that help to retrieve the content
of a *.xmf file:
   --> traverse_xmf
   --> point_xmf_name
   --> point_xmf_idx
   --> xmf_get_access

Other aux functions are:
   --> index_reverse
   --> subvol_process
   --> subvol_locate

#####################################################################################
Example of usage: Read full flow-field of a 3 block setup at timestep 500.

# your path to the FLOW files
path = "/home/username/.../sim_dir/FLOW/"

# Read blocks 1 to 3 in full.
grid = get_grids(wdir=path,block_ids=range(1,4))
flow = get_data(grid,wdir=path,tstep=500)

# Now we can access grid data to block 1 as grid[0].data of dimensions npx,npy,npz,3

# We can also access the flowfield of block 1 as flow[0].data of dimensions
# npx,npy,npz,5

# To access the grid and flow data from block 2 we would access grid[1].data and
# flow[1].data

#####################################################################################
"""

import warnings
import numpy as np
import struct as st
import xml.etree.ElementTree as ET
try:
    import lib_read_HiPSTAR_IO_shub as rf
    default_backend='Fortran'
except:
    ww = "lib_read_HiPSTAR could not find the read_plot3d module " \
                         +"(have you installed PLATUS?). File IO will be done in Python :("
    warnings.warn(ww)
    default_backend='Python'


class new_block(object):
    """
    Class to store single-block data.

    It also contains convenient methods to manipulate the data.
    """
    def __init__(self,blockid,arr,npx,npy,npz,idim,jdim,kdim,istart,iend,jstart,
                 jend,kstart,kend,iskip=1,jskip=1,kskip=1,subspace=False,
                 subspace_tag='',xdmf=None,timestep=None,time_span=None,
                 file_header=None,sampling_period=None,nstart=None,nend=None,
                 precision=None,backend=None,stat_bin=None,threed=True):

        self.blockid = blockid  # Block_id
        self.data = arr         # Data
        self.xdmf = xdmf        # xdmf pointer
        self.npx  = npx         # Global dataset dimensions
        self.npy  = npy
        self.npz  = npz
        self.idim  = idim       # Subset dimensions
        self.jdim  = jdim
        self.kdim  = kdim
        self.istart  = istart   # Global limits of the local subset
        self.iend    = iend
        self.jstart  = jstart
        self.jend    = jend
        self.kstart  = kstart
        self.kend    = kend
        self.iskip   = iskip
        self.jskip   = jskip
        self.kskip   = kskip
        self.subspace = subspace
        self.subspace_tag = subspace_tag
        self.stat_bin = stat_bin
        self.timestep = timestep
        self.time_span = time_span           # Time covered by the Sfile
        self.file_header = file_header
        self.sampling_period = sampling_period # Sampling period of the Sfile
        self.nstart          = nstart
        self.nend            = nend
        self.precision       = precision  # Reading precision
        self.backend         = backend    # Engine for reading (fortran of python)
        self.threed          = threed    # flag to determine 2D/3D simulation


        # Initialise lists to store Sfile convergence details. Only used if dealing with Sfiles
        # and the function averaged_Sfile_convergence.
        self.merged = [nstart] # List of Sfiles merged
        self.convergence = []
        for ii in range(arr.shape[-1]):
            self.convergence.append([np.nan])

        # Other variables
        self.purge_list = None

        # Aux - just in case you need to store stuff here.
        self.aux = None


    #########################
    ## The class functions ##
    #########################

    def update_file_header(self,_fh):
        """
        Updates the class variable file_header.
        """

        self.file_header = _fh

        return


    def update_period(self,_time_span):
        """
        Updates the time period after merging 2 Sfiles.
        """

        self.time_span += _time_span

        return


    def update_nend(self,_nend):
        """
        Updates the ending time step.
        """

        self.nend = _nend

        return


    def update_nsamples(self,_nsamples):
        """
        Updates the total number of samples after merging 2 Sfiles.
        """

        self.file_header[4] += _nsamples

        return


    def update_headers(self,Sfile2):
        """
        Updates the header and dt values after averaging 2 Sfiles.
        """

        self.update_period(Sfile2.time_span)
        self.update_nend(Sfile2.nend)
        self.update_nsamples(Sfile2.file_header[4])

        # Data capturing period (in time steps)
        ncapt = int(self.sampling_period/self.timestep)
        # Total number of time steps (not captured snapshots!)
        ntsteps = self.file_header[4]*ncapt

        # Update dt and sampling_step. They just need to be consistent with
        # the overall number of captured time steps and time period.
        self.timestep = float(self.time_span/ntsteps)
        self.sampling_period = float(self.timestep*ncapt)

        return


    def set_purge_list(self,_purge_list):
        """
        Sets the list of variables to be purged from the Sfiles.
        """

        self.purge_list = _purge_list

        return


    def mult_time(self):
        """
        Multiplies the time-averaged data by the sampling period.
        """

        self.data *= self.time_span

        return


    def div_time(self):
        """
        Divides the time-averaged data by the sampling period.
        """

        self.data /= self.time_span

        return


    def purge_header(self,ext_header):
        """
        Diffs the own file header (self.file_header) with an external one.
        After that, it works out the unmatched variables in ext_header.

        Then, it updates the file header, and produces a list of
        variables to be removed from the data array.

        This is of good use when averaging out 2 Sfiles with different variables.
        """

        _slice=[]
        fh1 = self.file_header
        fh2 = ext_header

        ngr = 0
        while ngr < fh1[5]:
            if fh1[6+ngr*2] not in fh2[6::2]:
                grp_offset = sum(fh1[7::2][:ngr])-1
                _slice.append(np.s_[grp_offset:grp_offset+fh1[7+ngr*2]])
                fh1[5] -= 1
                fh1 = np.delete(fh1,np.s_[6+ngr*2:8+ngr*2],axis=0)
                ngr = 0

            else:
                ngr += 1

        if fh1[5]==0:
            raise ValueError("Sfiles do not share any variable.")

        self.update_file_header(fh1)
        self.set_purge_list(_slice)

        return


    def purge_vars(self,Sfile2):
        """
        It purges the variables from self.data which do not match any variables
        of the object Sfile2. It ensures consistency when averaging multiple Sfiles
        which might have different variables stored.
        """

        # Check if headers are different. If not, just return.
        if not np.array_equal(self.file_header[5:],Sfile2.file_header[5:]):

            # Purge header and get the list of unmatched variables.
            self.purge_header(Sfile2.file_header)

            # Purge the data array of unmatched variables.
            for pg in self.purge_list:
                self.data = np.delete(self.data,pg,axis=3)


        return


    def process_Sfile(self,proc_ids="All",replace=True):

        """

        This function process the raw data from self.data and returs the post-processed
        flow statistics in the same array.

        The data is processed exaclty in the same way as it would be in HiPSTAR using the
        POSTPRO mode. From the information contained in self.file_header, we extract
        the groups of data contained in Sfile and process them accordingly. At the
        moment, only the processing of the following groups is available (processing
        IDs in brackets):

            - Group 0 [101] - Favre-averaged statistics (16 variables)
            - Group 1 [102] - Reynolds-averaged statistics (8 variables)
            - Group 4 [115] - LES statistics (1 variable)

        If required, feel free to implement more. If this is the case, bear in mind that
        the processing of the groups has to occur with increasing group ID order. For
        more info on the groups, have a look at the HiPSTAR wiki.

        The argument proc_ids defines which groups are processed. For example, if set to
        [102,115], only the statistics from groups 1 and 4 will be computed. The default
        is "All", which will return as many groups as available.

        *** Group IDs and processing IDs are not the same thing ***

        Notes: When computing the fluctuations of a quantity in a region where they are
               zero (e.g. velocity fluctuations at the wall), the square root might
               produce NaNs due to a negative value as an argument. This is purely due to
               numerical precision, and if the NaNs bother you, it is safe to put an
               np.abs inside the square root.

        """

        # Database with what's available for processing:
        # Format: [[groupid_0,nvars_in,[procid_0.0,nvars_out]],
        #          [groupid_1,nvars_in,[procid_1.0,nvars_out],
        #                              [procid_1.1,nvars_out]],
        #          ...
        #                  ]
        if self.threed:
              database = [
                  [0,16,[101,16]],  # Favre-averaged stats
                  [1,8 ,[102,8 ]],  # Reynolds-averaged stats
                  [4,1 ,[115,1 ]]   # Sub-grid viscosity
              ]
        else:
              database = [
                  [0,12,[101,12]],  # Favre-averaged stats
                  [1,6 ,[102,6 ]],  # Reynolds-averaged stats
                  [4,0 ,[115,0 ]]   # Sub-grid viscosity
              ]



        ## Fetch data from the Sfile.

        # Read file header and work out the groups in the Sfile.
        n_groups = self.file_header[5]
        if n_groups==0:
            raise ValueError("No groups available for processing.")
            return

        # Store group IDs and number of variables corresponding to each group
        grp_ids = np.zeros((n_groups,3),dtype=np.int)
        for i in range(n_groups):
            grp_ids[i,0] = self.file_header[5+(i*2+1)] # Group ID
            grp_ids[i,1] = self.file_header[5+(i*2+2)] # Number of variables in the group
            if i>0:
                grp_ids[i,2] = grp_ids[i-1,2]+grp_ids[i-1,1] # Varaible offset


        ## Check user options with what's in the Sfile, and also with what's implemented.

        # Setup output parameters.
        nvar = 0
        proc_out = []
        proc_off = []

        # Process all available groups.
        if proc_ids=="All":

            # Loop over the collected stat groups.
            for i in range(n_groups):
                # Loop over the available stat groups in the database.
                for group_avail in database:
                    # Look for a match in group ID.
                    if grp_ids[i,0]==group_avail[0]:
                        # Make sure the number of variables in the group are correct.
                        if grp_ids[i,1]!=group_avail[1]:
                            raise ValueError("Group %d does not have the expected number of variables." % (grp_ids[i,0]))
                        # If a match is found, do all the processing IDs from that stat group.
                        for procid in group_avail[2:]:
                            # Append processing ID
                            proc_out.append(procid[0])
                            # Append index offset
                            proc_off.append(grp_ids[i,2])
                            # Update total number of variables
                            nvar += procid[1]

                        break
                else:
                    print ("Processing of group with ID %d is not available." % (grp_ids[i,0]))


        # Process a subset of the availble groups.
        else:

            # Loop over the processing IDs
            for proc_id in proc_ids:
                 # Loop over the stat groups in the database.
                 for i,group_avail in enumerate(database):
                     # Loop over all the processing IDs from that stat group.
                     for procid in group_avail[2:]:
                         # Look for a match.
                         if proc_id==procid[0]:
                             # If the processing is available, now check if the stat group is in the Sfile.
                             for j in range(n_groups):
                                 if grp_ids[j,0]==group_avail[0]:
                                     break
                             else:
                                 raise ValueError("I can't do the processing with ID %d. Stat group %d is not in the Sfile." % (proc_id,group_avail[0]))

                             # Make sure the number of variables in the group are correct.
                             if grp_ids[j,1]!=group_avail[1]:
                                 raise ValueError("Group %d does not have the expected number of variables." % (grp_ids[i,0]))

                             # Append processing ID
                             proc_out.append(procid[0])
                             # Append index offset
                             proc_off.append(grp_ids[j,2])
                             # Update total number of variables
                             nvar += procid[1]

                             break
                     else:
                         continue
                     break
                 else:
                     print ("Processing ID %d is not available." % (proc_id))


        ## Done with the setup here. Allocate the arrays now and process the data.

        Sdata = self.data

        # Crate temp array
        (npx,npy,npz) = Sdata.shape[:3]
        Stemp = np.zeros((npx,npy,npz,nvar),dtype=np.float32)

        idx = 0
        # Loop over the processing IDs
        for i,ID in enumerate(proc_out):

            # Set offset - pointer to the beginning of the stat group.
            off = proc_off[i]

            ## Favre-averaged quantities
            if ID==101:

                Stemp[...,idx+0] =   Sdata[...,off+0]                                                      # 0  - Averaged density
                Stemp[...,idx+1] =   Sdata[...,off+2] / Sdata[...,off+0]                                   # 1  - Averaged u
                Stemp[...,idx+2] =   Sdata[...,off+3] / Sdata[...,off+0]                                   # 2  - Averaged v
                Stemp[...,idx+3] =   Sdata[...,off+4] / Sdata[...,off+0]                                   # 3  - Averaged T
                Stemp[...,idx+4] =   Sdata[...,off+6]                                                      # 4  - Averaged pressure
                Stemp[...,idx+5] =   Sdata[...,off+8]                                                      # 5  - Averaged molecular viscosity
                Stemp[...,idx+6] =   np.sqrt(Sdata[...,off+1]-Sdata[...,off+0]**2)                         # 6  - Density fluctuations
                Stemp[...,idx+7] =   np.sqrt(Sdata[...,off+5]/Sdata[...,off+0]-Stemp[...,off+3]**2)        # 7  - Temperature fluctuations
                Stemp[...,idx+8] =   np.sqrt(Sdata[...,off+7]-Sdata[...,off+6]**2)                         # 8  - Pressure fluctuations
                Stemp[...,idx+9] =  -Sdata[...,off+9]/Sdata[...,off+0]+Stemp[...,off+1]**2                 # 9  - tau_11 // u rms
                Stemp[...,idx+10] = -Sdata[...,off+11]/Sdata[...,off+0]+Stemp[...,off+1]*Stemp[...,off+2]  # 10 - tau_12 // uv rms
                Stemp[...,idx+11] = -Sdata[...,off+10]/Sdata[...,off+0]+Stemp[...,off+2]**2                # 11 - tau_22 // v rms
                if self.threed:
                    Stemp[...,idx+15] =  Sdata[...,off+15]/Sdata[...,off+0]                                    # 15 - Averaged w
                    Stemp[...,idx+12] = -Sdata[...,off+12]/Sdata[...,off+0]+Stemp[...,off+15]**2               # 12 - tau_33 // w rms
                    Stemp[...,idx+13] = -Sdata[...,off+13]/Sdata[...,off+0]+Stemp[...,off+1]*Stemp[...,off+15] # 13 - tau_13 // uw rms
                    Stemp[...,idx+14] = -Sdata[...,off+14]/Sdata[...,off+0]+Stemp[...,off+2]*Stemp[...,off+15] # 14 - tau_23 // vw rms
                    idx += 16 # Add the index offset
                else:
                    idx += 12 # Add the index offset


            ## Reynolds-averaged quantities
            elif ID==102:

                Stemp[...,idx]   =  Sdata[...,off]                                          # 12-16 - Averaged u
                Stemp[...,idx+1] =  Sdata[...,off+1]                                        # 13-17 - Averaged v
                Stemp[...,idx+2] =  Sdata[...,off+2]                                        # 14-18 - Averaged T
                Stemp[...,idx+3] =  np.sqrt(Sdata[...,off+3]-Sdata[...,off+2]**2)           # 15-19 - T rms
                Stemp[...,idx+4] =  np.sqrt(Sdata[...,off+4]-Sdata[...,off]**2)             # 16-20 - u rms
                Stemp[...,idx+5] =  np.sqrt(Sdata[...,off+5]-Sdata[...,off+1]**2)           # 17-21 - v rms
                if self.threed:
                    Stemp[...,idx+6] =  Sdata[...,off+6]                                        # 22 - Averaged w
                    Stemp[...,idx+7] =  np.sqrt(Sdata[...,off+7]-Sdata[...,off+6]**2)           # 23 - w rms
                    idx += 8
                else:
                    idx += 6

            elif ID==115:
                if self.threed:

                    ## LES
                    Stemp[...,idx] = Sdata[...,off]                                             # 24 - mu_sgs

                    idx += 1

            else:
                raise ValueError("The code should've never reached this point. There's a bug in the argument processing section.")

        if replace:
            del self.data
            self.data = Stemp
            return proc_out,nvar

        else:
            return Stemp



############################## AUXILIARY FUNCTIONS ##############################

def index_reverse(i,npp):
    """
    Create function to transform negative indices from the subvol argument in the
    functions below.
    """
    if i<0:
        i = npp+i
    return i

def bounds_check(subvol_bound,npt,iloc,idir,blid):

    if iloc==0:
       sloc = 'starting'
    elif iloc==1:
       sloc = 'ending'
    else:
       raise ValueError('sloc function not available.')

    if (np.abs(subvol_bound)>(npt-1)):
        raise ValueError('The Sub-volume specification for block %d has failed in ' % (blid) +
                         'the %s direction. The %s limit is out of bounds. ' % (idir,sloc))

    return

def subvol_process(subvol,block_ids):
    """
    Create a list (subvol_args) with subvol parameters for each block_id.
    """
    if subvol == None:
        subvol_args = [[None]]*len(block_ids)

    else:
        subvol_args = []
        subvol_match = False
        for i in range(len(block_ids)):
            for j in range(len(subvol)):
                for k in range(len(subvol[j][0])):
                    if block_ids[i]==subvol[j][0][k]:
                        subvol_args.append(subvol[j][1:])
                        subvol_match = True
                        break
                else:
                    continue
                break

            if not subvol_match:
                subvol_args.append([None])
            else:
                subvol_match = False

    return subvol_args


def subvol_locate(subvol_args_nb,npx,npy,npz,block_id):
    """
    Link subvolume options to global dimensions.
    """

    # Set defaults
    istart = 0
    iend   = npx-1
    iskip  = 1
    jstart = 0
    jend   = npy-1
    jskip  = 1
    kstart = 0
    kend   = npz-1
    kskip  = 1

    # Set the starting and ending indices, and also the point skipping.
    if subvol_args_nb != [None]:
        subvol_i = subvol_args_nb[0]
        subvol_j = subvol_args_nb[1]
        subvol_k = subvol_args_nb[2]

        if subvol_i != [None]:
            istart = index_reverse(subvol_i[0],npx)
            bounds_check(istart,npx,0,'x',block_id)
            if len(subvol_i)>1:
                iend = index_reverse(subvol_i[1],npx)
                bounds_check(iend,npx,1,'x',block_id)
                if len(subvol_i)>2:
                    iskip  = subvol_i[2]

        if subvol_j != [None]:
            jstart = index_reverse(subvol_j[0],npy)
            bounds_check(jstart,npy,0,'y',block_id)
            if len(subvol_j)>1:
                jend = index_reverse(subvol_j[1],npy)
                bounds_check(jend,npy,1,'y',block_id)
                if len(subvol_j)>2:
                    jskip = subvol_j[2]
        if subvol_k != [None]:
            kstart = index_reverse(subvol_k[0],npz)
            bounds_check(kstart,npz,0,'z',block_id)
            if len(subvol_k)>1:
                kend = index_reverse(subvol_k[1],npz)
                bounds_check(kend,npz,1,'z',block_id)
                if len(subvol_k)>2:
                    kskip = subvol_k[2]

    return istart,jstart,kstart,iend,jend,kend,iskip,jskip,kskip

################################# MAIN FUNCTIONS #################################

def get_grids(wdir='.',block_ids=[1],subvol=None,read_subspaces=False,
              subspaces_list=[],fmt='PLOT3D',precision='single',pointer=False,
              backend=default_backend,quiet=False):

    """
    This function returs the list of objects block_grid, which contains the grid
    coordinates from all blocks.

        - wdir: Path to the directory where the grid files are.

        - block_ids: List containing the ids of blocks to be read-in. If reading
                     consecutive blocks, call this function setting
                     nblocks=range(start,end). On the other hand, if reading
                     subspaces, this list is replaced by the data from the
                     subspaces_list argument.

                     THE IDs MUST MATCH THE ONES FROM HiPSTAR (i.e. first block has
                     block_id=1 and so on.)

        - subvol:    Specifies the start and end point in i,j and k direcitons for
                     the chosen blocks blocks. It is also possible to skip points
                     on each direction. See examples below.

                     The subvol argument is a list and its format is as follows:

            subvol = [
                      # 3D subvolume
                      [[blockid0],[istart,iend],[jstart,jend],[kstart,kend]],

                      # 3D subvolume skipping every jskip points along j
                      [[blockid1],[istart,iend],[jstart,jend,jskip],[kstart,kend]],

                      # Only cropping along j direction
                      [[blockid2],[None],[jstart,jend],[None]],

                      # Only cropping along i and j directions
                      [[blockid3],[istart,iend],[None],[kstart,kend]],

                      # Cropping along i and k for multiple blocks
                      [[blockid4,blockid5],[istart,iend],[jstart,jend],[None]],
                      ...
                     ]

        - read_subspaces: If True, this function will read the grid files from the
                          subspaces listed in subspaces_list. If False, function
                          will act as default and read grids from FLOW files.

        - subspaces_list: 2D list containing the names of the subspaces to be read-in,
                          and the block_id from the subspace. It supports multiple
                          block IDs per subspace.

        subspaces_list = [
                  # One subspace with name SStag0 located in block blockid0
                  [[blockid0],['SStag0']],

                  # Two subspaces with name SStag1 located in blocks blockid0 and blockid1
                                  [[blockid0,blockid1],['SStag1']],

                                  # Two subspaces with names SStag2 and SStag3 in block blockid2
                                  [[blockid2],['SStag2',SStag3]]
                 ]

        - fmt:       This argument defines the format of the grid files to be read
                     in. The available options are:

                         ** HiPSTAR: Reads the HiPSTAR grid files. This option is
                                     only available if read_subspaces=False. At the
                                     moment, this option only suports standard
                                     Cartesian grids (non -DGCYL).

                         ** PLOT3D:  Reads the *xyz PLOT3D files.

        - precision: Sets the reading precision. It defaults to 'single', but double
                     precision is also supported by setting precision='double'. This
                     option only affects PLOT3D format. It assumes that the file
                     header is with standard integers (not long ones!).

        - pointer:   If True, it will return the block_grid list whitout the actual
                     grid coordinates. This is useful in case we just want to read in
                     data with the function get_data and the grid coordinates are not
                     required.

        - backend:   It switches the reading routines for PLOT3D from 'Python' to
                     'Fortran'.

        - quiet: If True, it suppresses the print outputs.

    """


    # Check options are compatible.
    if read_subspaces:
        if len(subspaces_list) is 0:
            raise ValueError("subspaces_list is empty!")

        else:
            # SS block IDs & name tags
            block_ids = []
            subspaces_names = []
            SS=subspaces_list
            print(SS)
            for tag in SS[1]:
                for ID in SS[0]:
                    block_ids.append(ID)
                    subspaces_names.append(tag)
    else:
        subspaces_names = None

    # Define reading precision.
    if precision=='single':
        iprec=np.float32
        sprec='f'
        bprec=4
    elif precision=='double':
        iprec=np.float64
        sprec='d'
        bprec=8
    else:
        raise ValueError('Selected reading precision not available.')


    # Process subvol argument
    subvol_args = subvol_process(subvol,block_ids)

    block_grid = [] # Initialize the grid list

    # Block/Subspace loop
    for nb,block_id in enumerate(block_ids):

        ### Open file - also check if format is specified OK ###
        if not read_subspaces:

            if fmt=='HiPSTAR':
                # Read Blocks in HiPSTAR format.
                f1 = open(wdir+'z_grid_%d.dat' % (block_id))
                f2 = open(wdir+'r_grid_%d.dat' % (block_id))
                f3 = open(wdir+'span_grid.dat')

            elif fmt=='PLOT3D':
                # Read Blocks in PLOT3D
                thefile = wdir+'FLOW_phys_GRID_%d.xyz' % (block_id)
                f = open(thefile,'rb')

            else:
                raise ValueError('The format specified in argument the fmt is not supported for reading FLOW files.')

        else: # Subspaces

            if fmt=='PLOT3D':
                # Read Subspaces in PLOT3D
                thefile = wdir+subspaces_names[nb]+'_GRID_%d.xyz' % (block_id)
                f = open(thefile,'rb')

            else:
                raise ValueError('The format specified in argument the fmt is not supported for reading subspaces.')


        # Get the global extent of the block
        if fmt=='HiPSTAR':
            npx = int(f1.readline())
            npy = int(f2.readline())
            npz = int(f3.readline())

        else:
            # PLOT3D
            (npx,npy,npz,) = st.unpack('3i',f.read(3*4))


        # Now that we know the dimensions, check the subvolume options in i, j and k
        # directions.
        istart,jstart,kstart, \
        iend,jend,kend,       \
        iskip,jskip,kskip    = subvol_locate(subvol_args[nb],npx,npy,npz,block_id)

        idim  = len(range(istart,iend+1,iskip))
        jdim  = len(range(jstart,jend+1,jskip))
        kdim  = len(range(kstart,kend+1,kskip))
        npxy  = npx*npy
        npxyz = npxy*npz
        stride = bprec*(iskip-1)

        if (not quiet):
            print ("Global dimensions: %d, %d, %d ==> %d" % (npx,npy,npz,npx*npy*npz))
            print ("Subset dimensions: %d, %d, %d ==> %d" % (idim,jdim,kdim,idim*jdim*kdim))


        if not pointer:

            # Allocate grid
            grid = np.zeros((idim,jdim,kdim,3),dtype=iprec,order='FORTRAN')


            # Loop throught the coordinates and read the data.
            # We've got 2 readers here. First is the HiPSTAR one, followed by the PLOT3D
            # reader.

            if fmt=="HiPSTAR":

                ## Go to istart
                for i in range(istart):
                    f1.readline()

                ii = 0 ## Read x grid
                for i in range(istart,iend+1,iskip):
                    grid[ii,:,:,0] = float(f1.readline())

                    for skp in range(iskip-1): # Skipping
                       f1.readline()

                    ii += 1

                ## Go to jstart
                for i in range(jstart):
                    f2.readline()

                jj = 0 ## Read y grid
                for j in range(jstart,jend+1,jskip):
                    grid[:,jj,:,1] = float(f2.readline())

                    for skp in range(jskip-1): # Skipping
                       f2.readline()

                    jj += 1

                ## Go to kstart
                for i in range(kstart):
                    f3.readline()

                kk = 0 ## Read z grid
                for k in range(kstart,kend+1,kskip):
                    grid[...,kk,2] = float(f3.readline())

                    for skp in range(kskip-1): # Skipping
                       f3.readline()

                    kk += 1


            else:

                if backend=='Fortran':

                   # The file
                   iofile = thefile
                   # List of variables
                   fvars = np.array(range(3),dtype=int,order='FORTRAN')
                   # Global dimensions
                   abs_dims = np.array([npx,npy,npz],dtype=int,order='FORTRAN')
                   # Subvol argument
                   bounds = np.array([istart,iend,iskip,jstart,jend,jskip,kstart,kend,kskip]
                                          ,dtype=int,order='FORTRAN')
                   # Heather length
                   header_length = 12

                   if precision=='single':
                       rf.read_plot3d_single(iofile,grid,fvars,idim,jdim,kdim,3,
                                             abs_dims,bounds,header_length)

                   elif precision=='double':
                       rf.read_plot3d_double(iofile,grid,fvars,idim,jdim,kdim,3,
                                             abs_dims,bounds,header_length)


                elif backend=='Python':

                    for n in range(3):
                        kk = 0
                        for k in range(kstart,kend+1,kskip):
                            fbin = f.read(0)
                            for j in range(jstart,jend+1,jskip):
                                f.seek(12+(npxyz*n+npxy*k+npx*j+istart)*bprec)
                                for i in range(idim):
                                    fbin += f.read(bprec)
                                    f.seek(stride,1)

                            # Unpack the data
                            grid[:,:,kk,n] = np.reshape(st.unpack(('%d'+sprec)%(idim*jdim),fbin),
                                                       (idim,jdim),order='F')
                            kk += 1
                else:
                    raise ValueError('Unknown backend option.')


        else:
            # If running function in pointer mode, the grid array is empty.
            grid = None


        # Store grid in block_grid object list
        if read_subspaces:
            subspacetag = subspaces_names[nb]
        else:
            subspacetag = ""

        block_grid.append(new_block(block_id,grid,
                                    npx,npy,npz,
                                    idim,jdim,kdim,
                                    istart,iend,
                                    jstart,jend,
                                    kstart,kend,
                                    iskip,jskip,kskip,
                                    subspace = read_subspaces,
                                    subspace_tag = subspacetag,
                                    precision = precision,
                                    backend = backend))
        del grid

        # Close the files from the current block before moving on to the next one.
        if fmt=='HiPSTAR':
            f1.close()
            f2.close()
            f3.close()

        else:
            f.close()


    return block_grid


def get_data(block_grid,wdir='.',tstep=0,var=[],ss_var=0):

    """
    This function returs the list of objects block_data, which contains the required
    flow data from all blocks.

        - block_grid: Multiblock object list which contains the essential
                      information to read in the data. This argument MUST have been
                      generated with the function get_grid.

        - wdir: Path to the directory where the files are.

        - tstep: Time-step to be read in. THIS IS JUST A SINGLE INTEGER VALUE!

        - var: List specifying a subset of the variables available in the files. If
               empty, the var argument defaults to range(5) if reading flow files,
               or range(ss_var) if dealing with subspaces.

        - ss_var: If reading data from subspaces, ss_var is the total number of
                  variables stored in the subspace. When reading FLOW files, this
                  variable is not used.

    """

    # Check that options are correct.
    subspace = block_grid[0].subspace
    if subspace and ss_var==0:
        raise ValueError("If reading Subspaces, need to specify the total number of variables in the subspace!")

    if not subspace:
        if len(var)==0:
            var = range(5)
        if len(var)>5:
            raise ValueError("var is larger than the number of variables in the FLOW files.")

    else:
        if len(var)==0:
            var = range(ss_var)
        if len(var)>ss_var:
            raise ValueError("var is larger than the number of variables in the subspace files.")

    # Recover the block_ids list
    block_ids = []
    for i in range(len(block_grid)):
        block_ids.append(block_grid[i].blockid)


    block_data = [] # Initialise the block_data list

    for nb,block_id in enumerate(block_ids):

        # Get precision.
        precision = block_grid[nb].precision
        if precision=='single':
            iprec=np.float32
            sprec='f'
            bprec=4
        elif precision=='double':
            iprec=np.float64
            sprec='d'
            bprec=8
        else:
            raise ValueError('Selected reading precision not available.')

        # Get reading backend
        backend = block_grid[nb].backend


        # Pull the sub-volume pointers from the block_grid object.
        npx    = block_grid[nb].npx
        npy    = block_grid[nb].npy
        npz    = block_grid[nb].npz
        istart = block_grid[nb].istart
        iend   = block_grid[nb].iend
        jstart = block_grid[nb].jstart
        jend   = block_grid[nb].jend
        kstart = block_grid[nb].kstart
        kend   = block_grid[nb].kend
        iskip  = block_grid[nb].iskip
        jskip   = block_grid[nb].jskip
        kskip   = block_grid[nb].kskip


        idim  = len(range(istart,iend+1,iskip))
        jdim  = len(range(jstart,jend+1,jskip))
        kdim  = len(range(kstart,kend+1,kskip))
        npxy  = npx*npy
        npxyz = npxy*npz
        stride = bprec*(iskip-1)

        # Allocate the data
        data = np.zeros((idim,jdim,kdim,len(var)),dtype=iprec,order='FORTRAN')

        if block_grid[nb].subspace:
            thefile = wdir+block_grid[nb].subspace_tag+"_%d_var_%d_%d.raw" % (block_id,ss_var,tstep)
            f = open(thefile,'rb')
        else:
            thefile = wdir+"FLOW_phys_%d_%d.raw" % (block_id,tstep)
            f = open(thefile,'rb')


        # Read the data
        if backend=='Fortran':

            # The file
            iofile = thefile
            # List of variables
            fvars = np.array(var,dtype=int,order='FORTRAN')
            # Global dimensions
            abs_dims = np.array([npx,npy,npz],dtype=int,order='FORTRAN')
            # Subvol argument
            bounds = np.array([istart,iend,iskip,jstart,jend,jskip,kstart,kend,kskip]
                                   ,dtype=int,order='FORTRAN')
            # Heather length
            header_length = 28

            if precision=='single':
                rf.read_plot3d_single(iofile,data,fvars,idim,jdim,kdim,len(var),
                                      abs_dims,bounds,header_length)

            elif precision=='double':
                rf.read_plot3d_double(iofile,data,fvars,idim,jdim,kdim,len(var),
                                      abs_dims,bounds,header_length)


        elif backend=='Python':
            nidx = 0
            for n in var:
                kk = 0
                for k in range(kstart,kend+1,kskip):
                    fbin = f.read(0)
                    for j in range(jstart,jend+1,jskip):
                        f.seek(28+(npxyz*n+npxy*k+npx*j+istart)*bprec)
                        for i in range(idim):
                            fbin += f.read(bprec)
                            f.seek(stride,1)

                    data[:,:,kk,nidx] = np.reshape(st.unpack(('%d'+sprec)%(idim*jdim),fbin),
                                       (idim,jdim),order='F')
                    kk += 1
                nidx += 1

        else:
            raise ValueError('Unknown backend option.')


        # Store data in block_data object list
        block_data.append(new_block(block_id,data,
                                    npx,npy,npz,
                                    idim,jdim,kdim,
                                    istart,iend,
                                    jstart,jend,
                                    kstart,kend,
                                    iskip,jskip,kskip))
        del data
        f.close()

    return block_data



def getAndWrite_grids(wdir='.',wwdir='.',block_ids=[1],subvol=None,read_subspaces=False,
              subspaces_list=[],fmt='PLOT3D',precision='single',pointer=False,
              backend=default_backend,quiet=False):

    """
    This function returs the list of objects block_grid, which contains the grid
    coordinates from all blocks.

        - wdir: Path to the directory where the grid files are.

        - block_ids: List containing the ids of blocks to be read-in. If reading
                     consecutive blocks, call this function setting
                     nblocks=range(start,end). On the other hand, if reading
                     subspaces, this list is replaced by the data from the
                     subspaces_list argument.

                     THE IDs MUST MATCH THE ONES FROM HiPSTAR (i.e. first block has
                     block_id=1 and so on.)

        - subvol:    Specifies the start and end point in i,j and k direcitons for
                     the chosen blocks blocks. It is also possible to skip points
                     on each direction. See examples below.

                     The subvol argument is a list and its format is as follows:

            subvol = [
                      # 3D subvolume
                      [[blockid0],[istart,iend],[jstart,jend],[kstart,kend]],

                      # 3D subvolume skipping every jskip points along j
                      [[blockid1],[istart,iend],[jstart,jend,jskip],[kstart,kend]],

                      # Only cropping along j direction
                      [[blockid2],[None],[jstart,jend],[None]],

                      # Only cropping along i and j directions
                      [[blockid3],[istart,iend],[None],[kstart,kend]],

                      # Cropping along i and k for multiple blocks
                      [[blockid4,blockid5],[istart,iend],[jstart,jend],[None]],
                      ...
                     ]

        - read_subspaces: If True, this function will read the grid files from the
                          subspaces listed in subspaces_list. If False, function
                          will act as default and read grids from FLOW files.

        - subspaces_list: 2D list containing the names of the subspaces to be read-in,
                          and the block_id from the subspace. It supports multiple
                          block IDs per subspace.

        subspaces_list = [
                  # One subspace with name SStag0 located in block blockid0
                  [[blockid0],['SStag0']],

                  # Two subspaces with name SStag1 located in blocks blockid0 and blockid1
                                  [[blockid0,blockid1],['SStag1']],

                                  # Two subspaces with names SStag2 and SStag3 in block blockid2
                                  [[blockid2],['SStag2',SStag3]]
                 ]

        - fmt:       This argument defines the format of the grid files to be read
                     in. The available options are:

                         ** HiPSTAR: Reads the HiPSTAR grid files. This option is
                                     only available if read_subspaces=False. At the
                                     moment, this option only suports standard
                                     Cartesian grids (non -DGCYL).

                         ** PLOT3D:  Reads the *xyz PLOT3D files.

        - precision: Sets the reading precision. It defaults to 'single', but double
                     precision is also supported by setting precision='double'. This
                     option only affects PLOT3D format. It assumes that the file
                     header is with standard integers (not long ones!).

        - pointer:   If True, it will return the block_grid list whitout the actual
                     grid coordinates. This is useful in case we just want to read in
                     data with the function get_data and the grid coordinates are not
                     required.

        - backend:   It switches the reading routines for PLOT3D from 'Python' to
                     'Fortran'.

        - quiet: If True, it suppresses the print outputs.

    """


    # Check options are compatible.
    if read_subspaces:
        if len(subspaces_list) is 0:
            raise ValueError("subspaces_list is empty!")

        else:
            # SS block IDs & name tags
            print('boom')
            print(block_ids)
            print(subspaces_list)
            block_ids = []
            subspaces_names = []
            SS=subspaces_list
            print(SS)
            for tag in SS[1]:
                for ID in SS[0]:
                    block_ids.append(ID)
                    subspaces_names.append(tag)
    else:
        subspaces_names = None


    # Define reading precision.
    if precision=='single':
        iprec=np.float32
        sprec='f'
        bprec=4
    elif precision=='double':
        iprec=np.float64
        sprec='d'
        bprec=8
    else:
        raise ValueError('Selected reading precision not available.')


    # Process subvol argument
    subvol_args = subvol_process(subvol,block_ids)

    block_grid = [] # Initialize the grid list

    # Block/Subspace loop
    for nb,block_id in enumerate(block_ids):

        ### Open file - also check if format is specified OK ###
        if not read_subspaces:

            if fmt=='HiPSTAR':
                # Read Blocks in HiPSTAR format.
                f1 = open(wdir+'z_grid_%d.dat' % (block_id))
                f2 = open(wdir+'r_grid_%d.dat' % (block_id))
                f3 = open(wdir+'span_grid.dat')

            elif fmt=='PLOT3D':
                # Read Blocks in PLOT3D
                thefile = wdir+'FLOW_phys_GRID_%d.xyz' % (block_id)
                f = open(thefile,'rb')

            else:
                raise ValueError('The format specified in argument the fmt is not supported for reading FLOW files.')

        else: # Subspaces

            if fmt=='PLOT3D':
                # Read Subspaces in PLOT3D
                thefile = wdir+subspaces_names[nb]+'_GRID_%d.xyz' % (block_id)
#                wwdir  =  "/home/student.unimelb.edu.au/shubham1/Desktop/mattieu_data/LES_M05/WSubspaces/"
                wfile = wwdir+subspaces_names[nb]+'_GRID_%d.xyz' % (block_id)
                f = open(thefile,'rb')

            else:
                raise ValueError('The format specified in argument the fmt is not supported for reading subspaces.')


        # Get the global extent of the block
        if fmt=='HiPSTAR':
            npx = int(f1.readline())
            npy = int(f2.readline())
            npz = int(f3.readline())

        else:
            # PLOT3D
            (npx,npy,npz,) = st.unpack('3i',f.read(3*4))

        # Now that we know the dimensions, check the subvolume options in i, j and k
        # directions.
        istart,jstart,kstart, \
        iend,jend,kend,       \
        iskip,jskip,kskip    = subvol_locate(subvol_args[nb],npx,npy,npz,block_id)

        idim  = len(range(istart,iend+1,iskip))
        jdim  = len(range(jstart,jend+1,jskip))
        kdim  = len(range(kstart,kend+1,kskip))
        npxy  = npx*npy
        npxyz = npxy*npz
        stride = bprec*(iskip-1)

        if (not quiet):
            print ("Global dimensions: %d, %d, %d ==> %d" % (npx,npy,npz,npx*npy*npz))
            print ("Subset dimensions: %d, %d, %d ==> %d" % (idim,jdim,kdim,idim*jdim*kdim))


        if not pointer:

            # Allocate grid
            grid = np.zeros((idim,jdim,kdim,3),dtype=iprec,order='FORTRAN')


            # Loop throught the coordinates and read the data.
            # We've got 2 readers here. First is the HiPSTAR one, followed by the PLOT3D
            # reader.

            if fmt=="HiPSTAR":

                ## Go to istart
                for i in range(istart):
                    f1.readline()

                ii = 0 ## Read x grid
                for i in range(istart,iend+1,iskip):
                    grid[ii,:,:,0] = float(f1.readline())

                    for skp in range(iskip-1): # Skipping
                       f1.readline()

                    ii += 1

                ## Go to jstart
                for i in range(jstart):
                    f2.readline()

                jj = 0 ## Read y grid
                for j in range(jstart,jend+1,jskip):
                    grid[:,jj,:,1] = float(f2.readline())

                    for skp in range(jskip-1): # Skipping
                       f2.readline()

                    jj += 1

                ## Go to kstart
                for i in range(kstart):
                    f3.readline()

                kk = 0 ## Read z grid
                for k in range(kstart,kend+1,kskip):
                    grid[...,kk,2] = float(f3.readline())

                    for skp in range(kskip-1): # Skipping
                       f3.readline()

                    kk += 1


            else:

                if backend=='Fortran':

                   # The file
                   iofile = thefile
                   # List of variables
                   fvars = np.array(range(3),dtype=int,order='FORTRAN')
                   # Global dimensions
                   abs_dims = np.array([npx,npy,npz],dtype=int,order='FORTRAN')
                   # Subvol argument
                   bounds = np.array([istart,iend,iskip,jstart,jend,jskip,kstart,kend,kskip]
                                          ,dtype=int,order='FORTRAN')
                   # Heather length
                   header_length = 12

                   if precision=='single':  #  readwrite is only written for 'single' in fortran
                       rf.readwrite_plot3d_single(iofile,wfile,grid,fvars,idim,jdim,kdim,3,
                                             abs_dims,bounds,header_length)

                   elif precision=='double':
                       rf.read_plot3d_double(iofile,wfile,grid,fvars,idim,jdim,kdim,3,
                                             abs_dims,bounds,header_length)


                elif backend=='Python':
                    print('haha')

                    for n in range(3):
                        kk = 0
                        for k in range(kstart,kend+1,kskip):
                            fbin = f.read(0)
                            for j in range(jstart,jend+1,jskip):
                                f.seek(12+(npxyz*n+npxy*k+npx*j+istart)*bprec)
                                for i in range(idim):
                                    fbin += f.read(bprec)
                                    f.seek(stride,1)

                            # Unpack the data
                            grid[:,:,kk,n] = np.reshape(st.unpack(('%d'+sprec)%(idim*jdim),fbin),
                                                       (idim,jdim),order='F')
                            kk += 1
                else:
                    raise ValueError('Unknown backend option.')


        else:
            # If running function in pointer mode, the grid array is empty.
            grid = None


        # Store grid in block_grid object list
        if read_subspaces:
            subspacetag = subspaces_names[nb]
        else:
            subspacetag = ""

        block_grid.append(new_block(block_id,grid,
                                    npx,npy,npz,
                                    idim,jdim,kdim,
                                    istart,iend,
                                    jstart,jend,
                                    kstart,kend,
                                    iskip,jskip,kskip,
                                    subspace = read_subspaces,
                                    subspace_tag = subspacetag,
                                    precision = precision,
                                    backend = backend))
        del grid

        # Close the files from the current block before moving on to the next one.
        if fmt=='HiPSTAR':
            f1.close()
            f2.close()
            f3.close()

        else:
            f.close()


    return block_grid


def getAndWrite_data(block_grid,wdir='.',wwdir='.',tstep=0,var=[],ss_var=0):

    """
    This function returs the list of objects block_data, which contains the required
    flow data from all blocks.

        - block_grid: Multiblock object list which contains the essential
                      information to read in the data. This argument MUST have been
                      generated with the function get_grid.

        - wdir: Path to the directory where the files are.

        - tstep: Time-step to be read in. THIS IS JUST A SINGLE INTEGER VALUE!

        - var: List specifying a subset of the variables available in the files. If
               empty, the var argument defaults to range(5) if reading flow files,
               or range(ss_var) if dealing with subspaces.

        - ss_var: If reading data from subspaces, ss_var is the total number of
                  variables stored in the subspace. When reading FLOW files, this
                  variable is not used.

    """

    # Check that options are correct.
    subspace = block_grid[0].subspace
    if subspace and ss_var==0:
        raise ValueError("If reading Subspaces, need to specify the total number of variables in the subspace!")

    if not subspace:
        if len(var)==0:
            var = range(5)
        if len(var)>5:
            raise ValueError("var is larger than the number of variables in the FLOW files.")

    else:
        if len(var)==0:
            var = range(ss_var)
        if len(var)>ss_var:
            raise ValueError("var is larger than the number of variables in the subspace files.")

    # Recover the block_ids list
    block_ids = []
    for i in range(len(block_grid)):
        block_ids.append(block_grid[i].blockid)


    block_data = [] # Initialise the block_data list

    for nb,block_id in enumerate(block_ids):

        # Get precision.
        precision = block_grid[nb].precision
        if precision=='single':
            iprec=np.float32
            sprec='f'
            bprec=4
        elif precision=='double':
            iprec=np.float64
            sprec='d'
            bprec=8
        else:
            raise ValueError('Selected reading precision not available.')

        # Get reading backend
        backend = block_grid[nb].backend


        # Pull the sub-volume pointers from the block_grid object.
        npx    = block_grid[nb].npx
        npy    = block_grid[nb].npy
        npz    = block_grid[nb].npz
        istart = block_grid[nb].istart
        iend   = block_grid[nb].iend
        jstart = block_grid[nb].jstart
        jend   = block_grid[nb].jend
        kstart = block_grid[nb].kstart
        kend   = block_grid[nb].kend
        iskip  = block_grid[nb].iskip
        jskip   = block_grid[nb].jskip
        kskip   = block_grid[nb].kskip


        idim  = len(range(istart,iend+1,iskip))
        jdim  = len(range(jstart,jend+1,jskip))
        kdim  = len(range(kstart,kend+1,kskip))
        npxy  = npx*npy
        npxyz = npxy*npz
        stride = bprec*(iskip-1)

        # Allocate the data
        data = np.zeros((idim,jdim,kdim,len(var)),dtype=iprec,order='FORTRAN')

        if block_grid[nb].subspace:
            thefile = wdir+block_grid[nb].subspace_tag+"_%d_var_%d_%d.raw" % (block_id,ss_var,tstep)
#            wwdir  =  "/home/student.unimelb.edu.au/shubham1/Desktop/mattieu_data/LES_M05/WSubspaces/"
            wfile = wwdir+block_grid[nb].subspace_tag+"_%d_var_%d_%d.raw" % (block_id,ss_var,tstep)
            f = open(thefile,'rb')
        else:
            thefile = wdir+"FLOW_phys_%d_%d.raw" % (block_id,tstep)
            f = open(thefile,'rb')


        # Read the data
        if backend=='Fortran':

            # The file
            iofile = thefile
            # List of variables
            fvars = np.array(var,dtype=int,order='FORTRAN')
            # Global dimensions
            abs_dims = np.array([npx,npy,npz],dtype=int,order='FORTRAN')
            # Subvol argument
            bounds = np.array([istart,iend,iskip,jstart,jend,jskip,kstart,kend,kskip]
                                   ,dtype=int,order='FORTRAN')
            # Heather length
            header_length = 28

            if precision=='single':
                rf.readwrite_plot3d_single(iofile,wfile,data,fvars,idim,jdim,kdim,len(var),
                                      abs_dims,bounds,header_length)

            elif precision=='double':
                rf.read_plot3d_double(iofile,wfile,data,fvars,idim,jdim,kdim,len(var),
                                      abs_dims,bounds,header_length)


        elif backend=='Python':
            nidx = 0
            for n in var:
                kk = 0
                for k in range(kstart,kend+1,kskip):
                    fbin = f.read(0)
                    for j in range(jstart,jend+1,jskip):
                        f.seek(28+(npxyz*n+npxy*k+npx*j+istart)*bprec)
                        for i in range(idim):
                            fbin += f.read(bprec)
                            f.seek(stride,1)

                    data[:,:,kk,nidx] = np.reshape(st.unpack(('%d'+sprec)%(idim*jdim),fbin),
                                       (idim,jdim),order='F')
                    kk += 1
                nidx += 1

        else:
            raise ValueError('Unknown backend option.')


        # Store data in block_data object list
        block_data.append(new_block(block_id,data,
                                    npx,npy,npz,
                                    idim,jdim,kdim,
                                    istart,iend,
                                    jstart,jend,
                                    kstart,kend,
                                    iskip,jskip,kskip))
        del data
        f.close()

    return block_data


def traverse_xmf(obj,level = 0):
    """
    # This function goes through all the levels of an xmf file recursively.
    """
    tabspace = ''
    for i in range(level):
        tabspace +='\t'

    for subobj in obj:
        print (tabspace+'\033[91m'+str(subobj.tag)+"\033[0m"+str(subobj.attrib))
        traverse_xmf(subobj,level=level + 1)

    return


def point_xmf_name(obj,addr='',target='', start=True):
    """
    This function returns a new object pointing to a subelement specified with
    addr = '/obj_name/obj1_name/obj2_name...'

    The target attribute is usefull to get to a particular element when there
    are multiple objects of the same type at the same level.

    Watch out for non-uniqueness!
    """

    # Clean up addr on start.
    if start:
        addr = addr.split('/')
        while addr[0]=='':
            addr = addr[1:]
        while addr[-1]=='':
            addr = addr[:-1]

    rec = len(addr)

    for subobj in obj.iter(addr[0]):
        if rec>1 or target=='':
            obj = subobj
        else:
            if subobj.get('Name')==target:
                obj = subobj
                break

        if rec>1:
            obj = point_xmf_name(obj,target=target,addr=addr[1:],start=False)
            break

    return obj


def point_xmf_idx(obj,rec=0,idx=[0]):
    """
    This function returns a new object pointing to a subelement specified with
    the recursion level and the index at each level.
    """
    i = 0
    for subobj in obj:
        if i == idx[0]:
            obj = subobj

            if rec>1:
                obj = point_xmf_idx(obj,rec=rec-1,idx=idx[1:])

            break

        i += 1

    return obj


def xmf_get_access(obj):
    """
    Reads access and dataset parameters for reading.
    """
    seek_idx = int(obj.get('Seek'))   # Get possition in the binary file.
    prec = int(obj.get('Precision'))  # Dataset's precision.
    ntype = obj.get('NumberType')     # Data type
    if   prec==4 and ntype=="Float":
        btype = 'f'
    elif prec==8 and ntype=="Float":
        btype = 'd'
    elif prec==4 and ntype=="Integer":
        btype = 'i'
    else:
        raise ValueError('Type not known.')

    return seek_idx,prec,btype


def get_grids_stat_xmf(wdir='.',block_ids=[1],PRINT=False,subvol=None,backend=default_backend
                       ,stat_xmf='STAT.xmf',stat_bin='STAT_cont.bin',quiet=False):
    """
    This function returs the list of objects block_grid, which contains the grid
    coordinates from all blocks.

        - block_ids: List containing the IDs of the blocks to be read-in. This list
                     should follow the indexing given in HiPSTAR.

        - PRINT: When set to True, it enables printing the variables stored in the
                 STAT_cont.bin file.

        - subvol: specifies the start and end point in i,j and k direcitons the
                  blocks. See get_grids function for a more detailed description of
                  the format.

        - backend: Similarly to get_grids function, this argument selects the data
                   engine for reading. Options are 'Fortran' or 'Python'.

        - quiet: If True, it suppresses the print outputs.

    """

    import re

    filedir = wdir+stat_xmf
    tree = ET.parse(filedir)
    root = tree.getroot()

    # Process subvol argument
    subvol_args = subvol_process(subvol,block_ids)

    block_grid = [] # Initialize the grid object

    ### Read Grids ###
    thefile = wdir+stat_bin
    f = open(thefile)

    # Block loop
    for nb,block_id in enumerate(block_ids):

        # Get Block pointer in the xdmf tree.
        xdmf_block = point_xmf_name(root,addr='/Domain/Grid/Grid/',target='Block %d' % (block_id))

        # Get topology object to retrieve dimensions
        topo = point_xmf_name(xdmf_block,addr='/Grid/Topology')
        np_str = re.split('    |   |  | ',topo.get('Dimensions'))
        dim = len(np_str)
        npz = int(np_str[0])
        npy = int(np_str[1])
        npx = int(np_str[2])


        # Now that we know the dimensions, check the subvolume options in i, j and k
        # directions.
        istart,jstart,kstart, \
        iend,jend,kend,       \
        iskip,jskip,kskip    = subvol_locate(subvol_args[nb],npx,npy,npz,block_id)


        idim  = len(range(istart,iend+1,iskip))
        jdim  = len(range(jstart,jend+1,jskip))
        kdim  = len(range(kstart,kend+1,kskip))
        npxy = npx*npy


        # Allocate temp array and get the xdmf pointer
        grid = np.zeros((idim,jdim,kdim,dim),dtype=np.float32,order='FORTRAN')
        geo = point_xmf_name(xdmf_block,addr='/Grid/Geometry')

        if (not quiet):
            print ("Global dimensions: %d, %d, %d ==> %d" % (npx,npy,npz,npx*npy*npz))
            print ("Subset dimensions: %d, %d, %d ==> %d" % (idim,jdim,kdim,idim*jdim*kdim))


        # Loop throught the coordinates
        for n in range(dim):
            xdmf_data = point_xmf_idx(geo,idx=[n]) # Access the right xdmf element.
            seek_idx,prec,btype = xmf_get_access(xdmf_data) # Get binary, offset and precision
            stride = prec*(iskip-1)

            if prec!=4:
                raise ValueError('get_grids_stat_xmf is not ready for double precision data.')

            if backend=='Fortran':

                # The file
                iofile = thefile
                # List of variables - This is a dummy argument.
                fvars = np.array([0],dtype=int,order='FORTRAN')
                # Global dimensions
                abs_dims = np.array([npx,npy,npz],dtype=int,order='FORTRAN')
                # Subvol argument
                bounds = np.array([istart,iend,iskip,jstart,jend,jskip,kstart,kend,kskip]
                                       ,dtype=int,order='FORTRAN')
                # Heather length - For the stat_cont.bin files, this is the actual pointer
                # to the right variable we want to read in.
                header_length = seek_idx

                rf.read_plot3d_single(iofile,grid[...,n],fvars,idim,jdim,kdim,1,
                                      abs_dims,bounds,header_length)


            elif backend=='Python':
                kk = 0
                for k in range(kstart,kend+1,kskip):
                    fbin = f.read(0)
                    for j in range(jstart,jend+1,jskip):
                        f.seek(seek_idx+(npxy*k+npx*j+istart)*prec)
                        for i in range(idim):
                            fbin += f.read(prec)
                            f.seek(stride,1)

                    grid[:,:,kk,n] = np.reshape(st.unpack(('%d'+btype)%(idim*jdim),fbin),
                                                (idim,jdim),order='F')

                    kk += 1

        # Store grid in block_grid object list
        block_grid.append(new_block(block_id,grid,
                                    npx,npy,npz,
                                    idim,jdim,kdim,
                                    istart,iend,
                                    jstart,jend,
                                    kstart,kend,
                                    iskip,jskip,kskip,
                                    xdmf = xdmf_block,
                                    backend = backend,
                                    stat_bin = stat_bin))
        del grid

    f.close() # Close STAT_cont.bin file

    ## Print the quantities available
    if PRINT:
        i = 0
        for child in point_xmf_name(block_grid[0].xdmf,addr='/Grid'):
            if child.tag=='Attribute':
                print (i,(str(child.attrib.get('Name'))[1:-1]))
                i+=1

    return block_grid


def get_data_stat_xmf(block_grid,wdir='.',var=[]):
    """
    This function reads the data specified in var from all the blocks and returns
    it into the block_data object list.

        - block_grid: Multiblock object list which contains the essential
                      information to read the dataset. This argument MUST have been
                      genereated with the function get_grid_stat_xmf.

        - wdir: Path to the directory where the STAT_cont.bin and STAT.xmf are.

        - var: List of variables to be read-in from STAT_cont.bin file.
    """

    # Check that options are correct.
    if len(var) == 0:
        raise ValueError("var argument cannot be empty! Specify the variables you want to read-in.")

    # Recover the block_ids list
    block_ids = []
    for i in range(len(block_grid)):
        block_ids.append(block_grid[i].blockid)

    block_data = [] # Initialize data object

    ### Read Data ###

    thefile = wdir+block_grid[0].stat_bin
    f = open(thefile)

    for nb,block_id in enumerate(block_ids):

        # Pull the sub-volume pointers from the block_grid object.
        npx    = block_grid[nb].npx
        npy    = block_grid[nb].npy
        npz    = block_grid[nb].npz
        istart = block_grid[nb].istart
        iend   = block_grid[nb].iend
        jstart = block_grid[nb].jstart
        jend   = block_grid[nb].jend
        kstart = block_grid[nb].kstart
        kend   = block_grid[nb].kend
        iskip  = block_grid[nb].iskip
        jskip   = block_grid[nb].jskip
        kskip   = block_grid[nb].kskip

        idim  = len(range(istart,iend+1,iskip))
        jdim  = len(range(jstart,jend+1,jskip))
        kdim  = len(range(kstart,kend+1,kskip))
        npxy = npx*npy

        # Get the backend reader
        backend = block_grid[nb].backend

        # Allocate the data
        data = np.zeros((idim,jdim,kdim,len(var)),dtype=np.float32,order='FORTRAN')

        nidx = 0
        for n in var:

            # Get the info from the xmf
            xdmf_data = point_xmf_name(block_grid[nb].xdmf,addr='/Grid')
            xdmf_data = point_xmf_idx(xdmf_data,rec=2,idx=[2+n,0]) # Go one level further in the tree.
            seek_idx,prec,btype = xmf_get_access(xdmf_data) # Get binary, offset and precision
            stride = prec*(iskip-1)

            if prec!=4:
                raise ValueError('get_data_stat_xmf is not ready for double precision data.')


            if backend=='Fortran':

                # The file
                iofile = thefile
                # List of variables - This is a dummy argument.
                fvars = np.array([0],dtype=int,order='FORTRAN')
                # Global dimensions
                abs_dims = np.array([npx,npy,npz],dtype=int,order='FORTRAN')
                # Subvol argument
                bounds = np.array([istart,iend,iskip,jstart,jend,jskip,kstart,kend,kskip]
                                       ,dtype=int,order='FORTRAN')
                # Heather length - For the stat_cont.bin files, this is the actual pointer
                # to the right variable we want to read in.
                header_length = seek_idx

                rf.read_plot3d_single(iofile,data[...,nidx],fvars,idim,jdim,kdim,1,
                                      abs_dims,bounds,header_length)


            elif backend=='Python':
                kk = 0
                for k in range(kstart,kend+1,kskip):
                    fbin=f.read(0)
                    for j in range(jstart,jend+1,jskip):
                        f.seek(seek_idx+(npxy*k+npx*j+istart)*prec)
                        for i in range(idim):
                            fbin += f.read(prec)
                            f.seek(stride,1)

                    data[:,:,kk,nidx] = np.reshape(st.unpack(('%d'+btype)%(idim*jdim),fbin),
                                                   (idim,jdim),order='F')

                    kk += 1

            nidx += 1

        # Store data in block_data object list
        block_data.append(new_block(block_id,data,
                                    npx,npy,npz,
                                    idim,jdim,kdim,
                                    istart,iend,
                                    jstart,jend,
                                    kstart,kend,
                                    iskip,jskip,kskip,
                                    xdmf = block_grid[nb].xdmf))
        del data

    f.close() # Close the STAT_cont.bin file

    return block_data


def read_Sfile(block_grid,wdir='',tstep=0,quiet=False,group=[],var=[],threed=True):

    """
    This function returs the object list block_data, which contains the data from the
    Sfiles for all blocks.

        - block_grid: Multiblock object list which contains the essential information
                      to read in the data.
          This argument MUST have been generated with the function get_grid.

        - wdir: Path to the directory where the Sfile is.

        - tstep: Time-step to be read in. THIS IS JUST A SINGLE INTEGER VALUE! In
                 order to read multiple Sfiles at once, see the function
                 averaged_Sfile (or averaged_Sfile_convergence for more complex setups).

        - quiet: If True, it suppresses the print outputs.


        Two options to not read all the variables stored in an Sfile:

        - group: List containing a subset of groups to read from the Sfile. This list
                 should contain the group IDs of the required groups.

        - var: List containing a subset of variables to read from the Sfile. This list
               should contain the variable number in the order they were stored in the
               Sfile.

    """

    import os
    import copy

    # Check argument consistency.
    if ((group!=[])and(var!=[])):
        raise ValueError("Combination of input options not supported.")


    # Recover the block_ids list
    block_ids = []
    for i in range(len(block_grid)):
        block_ids.append(block_grid[i].blockid)


    block_data = [] # Initialise the block_data object

    for nb,block_id in enumerate(block_ids):

        # Pull the sub-volume pointers from the block_grid object.
        npx    = block_grid[nb].npx
        npy    = block_grid[nb].npy
        npz    = block_grid[nb].npz
        istart = block_grid[nb].istart
        iend   = block_grid[nb].iend
        jstart = block_grid[nb].jstart
        jend   = block_grid[nb].jend
        kstart = block_grid[nb].kstart
        kend   = block_grid[nb].kend
        iskip  = block_grid[nb].iskip
        jskip   = block_grid[nb].jskip
        kskip   = block_grid[nb].kskip

        idim  = len(range(istart,iend+1,iskip))
        jdim  = len(range(jstart,jend+1,jskip))
        kdim  = len(range(kstart,kend+1,kskip))
        npxy  = npx*npy
        npxyz = npxy*npz

        # Get the backend read in option.
        backend = block_grid[nb].backend

        thefile = wdir+"Sfile_b%d_%d" % (block_id,tstep)

        # Get the number of variables from the file size
        filesize = os.path.getsize(thefile)
        nvar = int(float(filesize-16)/float(npxyz*8))
        if (not quiet):
            print ("Sfile_b%d_%d contains %d variables." % (block_id,tstep,nvar))
        # and work out the size of the file header.
        header_size = filesize - (nvar*8*npxyz+16)


        f = open(thefile,'rb')

        # Read file header
        (dt,sampling_step) = st.unpack(('2d'),f.read(16))
        header = st.unpack('%di' % (header_size/4),f.read(header_size))


        # Get group info from the header.
        n_groups = header[5]
        # Create temp array with group IDs and number of variables corresponding to each group
        grp_ids = np.zeros((n_groups,3),dtype=np.int)
        for i in range(n_groups):
            grp_ids[i,0] = header[5+(i*2+1)] # Group ID
            grp_ids[i,1] = header[5+(i*2+2)] # Number of variables in the group
            if i>0:
                grp_ids[i,2] = grp_ids[i-1,2]+grp_ids[i-1,1] # Varaible offset

        # Make sure the number of variables is consistent with the file size.
        # If the error below is raised, you're likely to be trying to read an Sfile with
        # a block_grid object generated with the grid files from another case. This will
        # make the estimation of nvar based on the file size incorrect.
        if (nvar !=  grp_ids[:,1].sum()):
            raise ValueError("... or maybe not! Is block_grid consistent with the Sfile you're "+
                             "trying to read? Check your wdir!")
            return


        # If reading a subset of the variables or groups, process the header and
        # and create the list of variables to read. No need to do this if just reading
        # individual variables as the list is already provided as the var argument.
        if(group!=[]):

            # Duplicate the header
            header_new = []
            for h in range(5):
                header_new.append(header[h])
            # Add the number of groups to be read in
            header_new.append(len(group))
            # Add the group info
            grp_index = []
            for grp in group:
                for i in range(n_groups):
                    if grp==grp_ids[i,0]:
                        header_new.append(grp_ids[i,0]) # Group ID
                        header_new.append(grp_ids[i,1]) # Number of variables in thre group
                        grp_index.append(i)
                        break
                else:
                    print("Group ID %d is not present in the Sfile." % (grp))

            # Overwrite the header
            header = header_new

            # Create the list of variables.
            var=[]
            for idx in grp_index:
                for i in range(grp_ids[idx,1]):
                    var.append(grp_ids[idx,2]+i)


        elif(var!=[]):
            # Re-write the header with 0 groups.
            header_new = []
            for h in range(5):
                header_new.append(header[h])
            header_new.append(0)

            header = header_new




        # Two options to read the data here.
        # This is simply because of the way data is stored on the Sfiles,
        # it's just more efficient to keep the 2 routines for reading these in full,
        # or just a subset of the variables/groups.

        # Allocate the data
        if(var!=[]):
            nvar_full = nvar
            nvar = len(var)
        data = np.zeros((idim,jdim,kdim,nvar),dtype=np.float64,order='FORTRAN')

        # Read the full thing.
        if (var==[]):
            if backend=='Fortran':
                # The file
                iofile = thefile
                # List of variables - With the Sfiles, we read all of them.
                fvars = np.array(range(nvar),dtype=int,order='F')
                # Global dimensions
                abs_dims = np.array([npx,npy,npz],dtype=int,order='F')
                # Subvol argument
                bounds = np.array([istart,iend,iskip,jstart,jend,jskip,kstart,kend,kskip]
                                       ,dtype=int,order='F')
                # Header length
                header_length = 16+header_size

                rf.read_sfile(iofile,data,fvars,idim,jdim,kdim,nvar,
                                  abs_dims,bounds,header_length)


            elif backend=='Python':
                kk = 0
                for k in range(kstart,kend+1,kskip):
                    jj = 0
                    for j in range(jstart,jend+1,jskip):
                        fbin = f.read(0)
                        for i in range(istart,iend+1,iskip):
                            f.seek(16+header_size+(npxy*k+npx*j+i)*8*nvar)
                            fbin += f.read(8*nvar)

                        data[:,jj,kk,:] = np.reshape(st.unpack(('%dd')%(nvar*idim),fbin),
                                                     (idim,nvar),order='C')

                        jj += 1
                    kk += 1


        # Read a subset of the variables.
        else:
            if backend=='Fortran':
                # The file
                iofile = thefile
                # List of variables - With the Sfiles, we read all of them.
                fvars = np.array(var,dtype=int,order='F')
                # Global dimensions
                abs_dims = np.array([npx,npy,npz],dtype=int,order='F')
                # Subvol argument
                bounds = np.array([istart,iend,iskip,jstart,jend,jskip,kstart,kend,kskip]
                                       ,dtype=int,order='F')
                # Header length
                header_length = 16+header_size

                rf.read_sfile_vars(iofile,data,fvars,idim,jdim,kdim,nvar,nvar_full,
                                  abs_dims,bounds,header_length)


            elif backend=='Python':
                kk = 0
                for k in range(kstart,kend+1,kskip):
                    jj = 0
                    for j in range(jstart,jend+1,jskip):
                        fbin = f.read(0)
                        for i in range(istart,iend+1,iskip):
                            for n in var:
                                f.seek(16+header_size+((npxy*k+npx*j+i)*nvar_full+n)*8)
                                fbin += f.read(8)

                        data[:,jj,kk,:] = np.reshape(st.unpack(('%dd')%(nvar*idim),fbin),
                                                     (idim,nvar),order='C')

                        jj += 1
                    kk += 1


        # Store data in block_data object list
        block_data.append(new_block(block_id,data,
                                    npx,npy,npz,
                                    idim,jdim,kdim,
                                    istart,iend,
                                    jstart,jend,
                                    kstart,kend,
                                    iskip,jskip,kskip,
                                    file_header=np.array(header),
                                    timestep=dt,
                                    sampling_period = sampling_step,
                                    time_span = sampling_step*header[4],
                                    nstart = tstep,
                                    nend = tstep,threed=threed))
        del data
        f.close()

    return block_data


def averaged_Sfile(block_grid,wdir,tsteps,Sfile=None,return_last=False,group=[],var=[],quiet=False,threed=True):

    """
    This function returs the object Sfile, which contains the averaged data from
    multiple Sfiles fror all blocks. There is the posibility to extend the average of
    an existing Sfile object (see sfile argument). In order to get further info on the
    convergence of the averaging, consider using averaged_Sfile_convergence instead.

        - block_grid: Multiblock grid object which contains the essential information
                      to read in the data. This argument MUST have been generated
                      with the function get_grid.

        - wdir: Path to the directory where the Sfiles are.

        - tsteps: List containing the timestep of the Sfiles to be added on the final
                  average.

        - Sfile: Sfile object obtained with either read_Sfile or averaged_Sfile
                 functions. If you read in an Sfile with any of these functions, you
                 can just extend the averaging period of this object with extra sfiles.

        - return_last: if true, it also returns the last Sfile_temp object.

        - group: see read_Sfile.

        - var: see read_Sfile.

        - quiet: see read_Sfile.

    """

    # Set defaults
    Sfile2 = None

    # Recover the number of blocks
    nblocks = len(block_grid)

    if (Sfile==None):
        # Read the first Sfile to allocate the memory.
        Sfile = read_Sfile(block_grid,wdir,tstep=tsteps[0],group=group,var=var,quiet=quiet,threed=threed)
        div_stat = Sfile[0].time_span
    else:
        # Prepend a dummy timestep if the first Sfile was passed as an argument.
        tsteps = [0]+tsteps

    for nb in range(nblocks):
        Sfile[nb].mult_time()

    if len(tsteps)>1:
        for nit in tsteps[1:]: # Now just add the contribution of the other files.

            Sfile_temp = read_Sfile(block_grid,wdir,tstep=nit,group=group,var=var,quiet=quiet,threed=threed)

            for nb in range(nblocks):
                # Purge unmatched variables - won't do anything if the file headers match.
                Sfile[nb].purge_vars(Sfile_temp[nb])
                Sfile_temp[nb].purge_vars(Sfile[nb])
                # Add the stat data together
                Sfile_temp[nb].mult_time()
                Sfile[nb].data += Sfile_temp[nb].data
                # Update headers
                Sfile[nb].update_headers(Sfile_temp[nb])
                # Update merged list
                Sfile[nb].merged.append(nit)

            # Do not delete the Sfile_temp if required on return
            if ((return_last==False) or (return_last and (nit!=tsteps[-1]))):
                del Sfile_temp
            else:
                # If last Sfile is required, recover the local average and set output pointer.
                for nb in range(nblocks):
                    Sfile_temp[nb].div_time()
                Sfile2 = Sfile_temp

    for nb in range(nblocks):
        Sfile[nb].div_time()

    if(Sfile2!=None):
        return Sfile,Sfile2
    else:
        return Sfile


def averaged_Sfile_convergence(block_grid,wdir,tsteps,Sfile=None,return_last=False,group=[],var=[]
                              ,err_processing=None,quiet=False,threed=True):

    """
    This function is similar to averaged_Sfile, but it has the additional feature that tracks the
    statistical convergence of a dataset on the fly as the Sfiles are merged together.

    The output: The list with the individual Sfiles that have been merged so far can be is obtained in
    Sfile[:].merged. Also the list Sfile[:].convergence contains the evolution of the maximum rate of
    change from all points for each quantity as the Sfiles are added together. This is computed as

                         convergence_of_var_{n} = (var_{n-1}-var_{n})/var_{n-1},

    where n is the Sfile index.


        - err_processing: User-defined function that process the Sfiles raw-data. This option is very
                          useful just in case someone reads a subset of variables from the Sfiles and
                          want to check the convergence of a derived quantity from the raw variables.
                          In such case, this function should compute such derived quantity. See the
                          dummy_function defined below to see the arguments/returns this function
                          requires. If no err_processing function is defined, the dummy_function is
                          used, which will simply track the convrgence of the quantities are stored
                          on the Sfiles.

        For all the other arguments, please refer to averaged_Sfile.

    """
    import copy

    # Dummy function
    def dummy_function(data):
        var_pointers = []
        for ii in range(data.shape[-1]):
            var_pointers.append(data[...,ii])
        return var_pointers

    # Create function pointer if none passed as argument.
    if err_processing==None:
        err_processing=dummy_function


    if ((Sfile==None) and (len(tsteps)<2)):
        raise ValueError("It doesn't make sense to call this function for just one Sfile.")
        return

    # Add dummy timestep to the list
    elif (len(tsteps)<2):
        tsteps = [0]+tsteps


    if (Sfile==None):
        Sfile = read_Sfile(block_grid,wdir,tstep=tsteps[0],group=group,var=var,threed=threed)

    # Get number of variables
    nvar = Sfile[0].data.shape[-1]

    for stp in tsteps[1:]:

        # Sfile_prev is the current status of the total average.
        Sfile_prev = copy.deepcopy(Sfile)
        # Sfile is the updated total average and Sfile_last is literally the last Sfile.
        Sfile,Sfile_last = averaged_Sfile(block_grid,wdir,tsteps=[stp],Sfile=Sfile,group=group,var=var
                                               ,return_last=True)

        # Create temporary error list
        block_error = []
        for nb in range(nvar):
            block_error.append([])

        # Append max difference per block & variable
        for nb in range(len(Sfile)):
            e1 = err_processing(Sfile_prev[nb].data)
            e2 = err_processing(Sfile[nb].data)
            for i in range(len(e1)):
                block_error[i].append(np.max(np.abs(e1[i]-e2[i])/e1[i]))
                
        # Get max error from all blocks and append to convergence
        # (All blocks contain the same info)
        for nb in range(len(Sfile)):
            for i,error in enumerate(Sfile[nb].convergence):
                error.append(max(block_error[i]))

    if return_last:
        return Sfile,Sfile_last
    else:
        return Sfile


def process_Sfile(Sfile,process_groups="All"):

    """

    Wrapper to process the raw data from Sfile, returning the post-processed flow
    statistics in the same structure Sfile.

        - Sfile: This argument contains the multi-block data from the HiPSTAR Sfiles.
                 See read_Sfile for more info on what this data structure contains.
                 This argument MUST have been generated with the functions read_Sfile
                 or averaged_Sfile.

        - process_groups: See process_Sfile function in the new_block class.

    See the function process_Sfile in the class new_block for more information.
    """

    for block in Sfile:
        post_ids,nv = block.process_Sfile(proc_ids=process_groups)

    n_post_groups = len(post_ids)
    out_statement = "Sfiles with groups "
    for i in range(n_post_groups):
        out_statement += "%d, " % (post_ids[i])
    out_statement += "were processed correctly. Total number of quantities: %d" % (nv)

    print (out_statement)

    return Sfile


def write_Sfile(Sfile,dirout,offset=0):

    """

    This function takes the data from the Sfile object list and writes it out as a
    HiPSTAR Sfile.

        - Sfile: This argument contains the Sfile data, which MUST have
                 been generated with the function read_Sfile or averaged_Sfile.

        - dirout: Path to the directory where the file is to be written.

        - offset: If set to 0, it will open a new file and write the headers,
                  alongside the first subdomain. If different than 0, it will just
                  append subsequent subdomains to the file.

    """


    for block in Sfile:

        thefile = dirout+"Sfile_b%d_%d" % (block.blockid,block.nend)


        # Write header
        if offset == 0:
            f = open(thefile,'wb')
            head_ = [block.timestep,block.sampling_period]
            f.write(st.pack('2d',*head_))
            f.write(st.pack('%di' % (len(block.file_header)),*block.file_header))

        else:
            f = open(thefile,'ab')

        (npx,npy,npz,nvar) = block.data.shape
        npt = npx*nvar

        for k in range(npz):
            for j in range(npy):
                f.write(st.pack('%dd' % (npt),*np.ndarray.flatten(block.data[:,j,k,:],order='C')))

        f.close()

    return


def merge_Sfile(grid_dir='',sfile_dir='',out_dir='',block_ids=[1],tsteps=[0],fmt='PLOT3D'):

    """
    This function merges multiple Sfiles into one. This has the benefit of using up
    less disk space. In particular, this function is designed with huge Sfiles in
    mind, where one could not use the averaged_Sfile function in the whole domain.
    Thus, here we extract x-y slices of the Sfiles and we average them and write
    them to the new file. The subsequent slices are appended to this merged Sfile.
    This is of course much slower than doing it all at once, so you're only advised
    to use this function when the whole Sfile does not fit in your computer's
    memory. Otherwise you can merge smaller Sfiles calling averaged_Sfile and then
    write the results out with write_Sfile.

        - grid_dir: Directory where the case grid files are. See fmt argument.

        - sfile_dir: Directory where the Sfiles are.

        - out_dir: Directory where the merged Sfile will be written to.

        - block_ids: List of block ids from which we want to average their Sfiles.
                     Bear in mind that these IDs are HiPSTAR IDs, being the first
                     block ID=1, second block ID=2 and so on.

        - tsteps: List containing the timestep of the Sfiles to be averaged.

        - fmt: Format of the grid files. See get_grids function for all the options.

    """


    # Get absoulte dimensions.
    abs_grids = get_grids(grid_dir,block_ids,fmt=fmt,pointer=True)

    # Assume that all blocks have the same npz
    for k in range(abs_grids[0].npz):

        print (">> Merging Sfiles... %d" % (int((float(k)/float(abs_grids[0].npz))*100))+"%")

        # Assemble the subvol argument - this defines how we split the file
        # merge into smaller chunks.
        sv=[]
        for block in abs_grids:
            sv.append([[block.blockid],None,None,[k,k]])

        # Get the subvolume grid pointers
        sub_grids = get_grids(grid_dir,block_ids,subvol=sv,fmt=fmt,pointer=True)

        # Get the averaged Sfile subvolume.
        sub_sfile = averaged_Sfile(sub_grids,sfile_dir,tsteps)

        # Write files
        write_Sfile(sub_sfile,out_dir,offset=k)

    print (">> Done!")

    return


def write_STAT_cont(Sfile,block_grid,dirout,cgns=False,stat_bin='STAT_cont.bin'):

    """

    This function takes the data from the Sfile object list and creates the
    corresponding STAT_cont.bin file.

        - Sfile: This argument contains the processed statistics, which MUST have
                 been generated with the function process_Sfile.

        - block_grid: Multiblock grid object which contains the grid coordinates of
                      the data contained in Sfile. This argument MUST have been
                      generated with the function get_grid.

        - dirout: Path to the directory where the file is to be written.

        - cgns: NOT DOING ANYTHING YET. In the future, this function might support
                output in cgns format. Feel free to extend it yourself :)

    The header of the STAT_cont.bin file contains information about the groups
    present in the file. Each group has its own particular ID, which is different to
    those specified earlier in the function process_Sfile. Currently only the the
    following groups are available:

        - Group 0 - Favre-averaged statistics (16 variables)
        - Group 1 - Reynolds-averaged statistics (8 variables)
        - Group 4 - LES statistics (1 variable)

    For more info on the groups, have a look at the HiPSTAR wiki.

    """


    nblocks = len(Sfile)

    # Read file header and work out the groups to be written out.
    n_post_groups = Sfile[0].file_header[5]

    # Store group IDs and number of variables corresponding to each group
    post_ids = np.zeros((n_post_groups,2),dtype=np.int)
    for i in range(n_post_groups):
        post_ids[i,0] = Sfile[0].file_header[5+(i*2+1)] # Group ID
        post_ids[i,1] = Sfile[0].file_header[5+(i*2+2)] # Number of variables in the group

    # Build file header
    Stat_var = post_ids[:,1].sum()
    n_head = 5+3*nblocks+2*(n_post_groups+1)+3

    ibuf = np.zeros((n_head),dtype=np.int32)

    # Stat_cont_bin Header
    if Sfile[0].threed:
        ibuf[0] = 1  # 0 = 2D - 1 = 3D
    else:
        ibuf[0] = 0  # 0 = 2D - 1 = 3D
    ibuf[1] = 1 # 1 - Single precision output
    ibuf[2] = 1 # 0 - Single block (not implem.) - 1 multi-block
    ibuf[3] = nblocks
    for nb in range(nblocks):
        ibuf[nb*3+4] = Sfile[nb].idim
        ibuf[nb*3+5] = Sfile[nb].jdim
        ibuf[nb*3+6] = Sfile[nb].kdim

    i = 4+nblocks*3
    ibuf[i] = n_post_groups+1
    ibuf[i+1] = 100
    ibuf[i+2] = 3

    i += 2
    for nn in range(n_post_groups):
        # Transform the HiPSTAR group IDs to the ones supported by the PLATUS library.
        if post_ids[nn,0]==0:
            idout=101
        elif post_ids[nn,0]==1:
            idout=102
        elif post_ids[nn,0]==4:
            idout=115
        else:
            raise ValueError("Group ID %d not known." % (post_ids[nn,0]))

        ibuf[i+1] = idout
        ibuf[i+2] = post_ids[nn,1]
        i += 2

    ibuf[i+1] = int(Sfile[0].time_span/Sfile[0].sampling_period) # number of samples
    ibuf[i+2] = Sfile[0].nstart
    ibuf[i+3] = Sfile[0].nend


    with open(dirout+stat_bin,'wb') as f:

        # Write header
        f.write(st.pack('%di' % (n_head),*ibuf[:]))
        f.write(st.pack('f',Sfile[0].time_span))

        for nb in range(nblocks):

            npt = Sfile[nb].idim*Sfile[nb].jdim

            # Write Grid
            for n in range(3):
                for k in range(Sfile[nb].kdim):
                    f.write(st.pack('%df' % (npt),*np.ndarray.flatten(np.float32(block_grid[nb].data[:,:,k,n]),order='F')))

            # Write Data
            for n in range(Stat_var):
                for k in range(Sfile[nb].kdim):
                    f.write(st.pack('%df' % (npt),*np.ndarray.flatten(np.float32(Sfile[nb].data[:,:,k,n]),order='F')))



def grid_metrics(grid_obj,ord_deriv=4):
    """
    This function computes the grid metric terms required to work with derivatives
    in generalised 3D cartesian coordinates. Derivatives are computed with the deriv
    library, which contains the same numerical schemes which are used in HiPSTAR.
    Default is 4th order accuracy (see deriv.py for more details).

    This function returns the multi-block metrics object.
    """

    import deriv
    dv = deriv.deriv()

    # Declare the object structure.
    metric_obj = []

    class three_d_metric_obj(object):
        def __init__(self,x_xi,y_xi,z_xi,x_et,y_et,z_et,x_th,y_th,z_th,
                     Jdet_xiet,Jdet_xith,Jdet_etth,Jdet):
            self.x_xi = x_xi
            self.y_xi = y_xi
            self.z_xi = z_xi
            self.x_et = x_et
            self.y_et = y_et
            self.z_et = z_et
            self.x_th = x_th
            self.y_th = y_th
            self.z_th = z_th
            self.Jdet_xiet = Jdet_xiet
            self.Jdet_xith = Jdet_xith
            self.Jdet_etth = Jdet_etth
            self.Jdet = Jdet

    # Loop through the blocks and compute the metrics for each one.
    nblocks = len(grid_obj)

    for nb in range(nblocks):

        # Work with the sub-volume dimensions.
        npx = grid_obj[nb].idim
        npy = grid_obj[nb].jdim
        npz = grid_obj[nb].kdim
        grid = grid_obj[nb].data

        x_xi = np.zeros((npx,npy,npz)) # 1D terms
        x_et = np.zeros((npx,npy,npz))
        x_th = np.zeros((npx,npy,npz))
        y_xi = np.zeros((npx,npy,npz))
        y_et = np.zeros((npx,npy,npz))
        y_th = np.zeros((npx,npy,npz))
        z_xi = np.zeros((npx,npy,npz))
        z_et = np.zeros((npx,npy,npz))
        z_th = np.zeros((npx,npy,npz))
        Jdet_xiet = np.zeros((npx,npy,npz)) # 2D Determinants
        Jdet_xith = np.zeros((npx,npy,npz))
        Jdet_etth = np.zeros((npx,npy,npz))
        Jdet = np.zeros((npx,npy,npz)) # 3D Determinant

        for j in range(npy):
            x_xi[:,j,:] = dv.d2_carpenter(grid[:,j,:,0],order=ord_deriv)
            y_xi[:,j,:] = dv.d2_carpenter(grid[:,j,:,1],order=ord_deriv)
            z_xi[:,j,:] = dv.d2_carpenter(grid[:,j,:,2],order=ord_deriv)
        for i in range(npx):
            x_et[i,:,:] = dv.d2_carpenter(grid[i,:,:,0],order=ord_deriv)
            y_et[i,:,:] = dv.d2_carpenter(grid[i,:,:,1],order=ord_deriv)
            z_et[i,:,:] = dv.d2_carpenter(grid[i,:,:,2],order=ord_deriv)
        for i in range(npx):
            x_th[i,:,:] = dv.d2_carpenter(grid[i,:,:,0],order=ord_deriv,axis=1)
            y_th[i,:,:] = dv.d2_carpenter(grid[i,:,:,1],order=ord_deriv,axis=1)
            z_th[i,:,:] = dv.d2_carpenter(grid[i,:,:,2],order=ord_deriv,axis=1)

        for k in range(npz):
            Jdet_xiet[...,k] = x_xi[...,k]*y_et[...,k] - x_et[...,k]*y_xi[...,k]
        for j in range(npy):
            Jdet_xith[:,j,:] = x_xi[:,j,:]*z_th[:,j,:] - x_th[:,j,:]*z_xi[:,j,:]
        for i in range(npx):
            Jdet_etth[i,...] = y_et[i,...]*z_th[i,...] - y_th[i,...]*z_et[i,...]

        Jdet = x_xi*y_et*z_th + x_et*y_th*z_xi + x_th*y_xi*z_et \
              -x_th*y_et*z_xi - x_et*y_xi*z_th - x_xi*y_th*z_et

        # Store the metrics.
        metric_obj.append(three_d_metric_obj(x_xi,y_xi,z_xi,x_et,y_et,z_et,x_th,y_th,z_th \
                                         ,Jdet_xiet,Jdet_xith,Jdet_etth,Jdet))


    return metric_obj

def Q(arr,grid,ord_deriv=4):

    import deriv as dv
    dv = dv.deriv(norm='diag')

    # Get the metrics
    if grid.aux == None:
        tdset=[]
        tdset.append(grid)
        grid.aux = grid_metrics(tdset,ord_deriv=ord_deriv)[0]

    dxdxi = grid.aux.x_xi
    dydxi = grid.aux.y_xi
    dydet = grid.aux.y_et
    dxdet = grid.aux.x_et
    dzdth = grid.aux.z_th
    Jdet_xiet = grid.aux.Jdet_xiet

    dudxi = np.zeros(arr.shape[:-1],dtype=np.float32)
    dvdxi = np.zeros(arr.shape[:-1],dtype=np.float32)
    dwdxi = np.zeros(arr.shape[:-1],dtype=np.float32)
    dudet = np.zeros(arr.shape[:-1],dtype=np.float32)
    dvdet = np.zeros(arr.shape[:-1],dtype=np.float32)
    dwdet = np.zeros(arr.shape[:-1],dtype=np.float32)
    dudz = np.zeros(arr.shape[:-1],dtype=np.float32)
    dvdz = np.zeros(arr.shape[:-1],dtype=np.float32)
    dwdz = np.zeros(arr.shape[:-1],dtype=np.float32)

    (npx,npy,npz) = arr.shape[:-1]

    for j in range(npy):
        dudxi[:,j,:] = dv.d2_carpenter(arr[:,j,:,1],order=ord_deriv)
        dvdxi[:,j,:] = dv.d2_carpenter(arr[:,j,:,2],order=ord_deriv)
        dwdxi[:,j,:] = dv.d2_carpenter(arr[:,j,:,3],order=ord_deriv)
    for i in range(npx):
        dudet[i,:,:] = dv.d2_carpenter(arr[i,:,:,1],order=ord_deriv)
        dvdet[i,:,:] = dv.d2_carpenter(arr[i,:,:,2],order=ord_deriv)
        dwdet[i,:,:] = dv.d2_carpenter(arr[i,:,:,3],order=ord_deriv)
    for i in range(npx):
        dudz[i,:,:] = dv.d2_carpenter(arr[i,:,:,1],order=ord_deriv,axis=1)*(1.0/dzdth[i,:,:])
        dvdz[i,:,:] = dv.d2_carpenter(arr[i,:,:,2],order=ord_deriv,axis=1)*(1.0/dzdth[i,:,:])
        dwdz[i,:,:] = dv.d2_carpenter(arr[i,:,:,3],order=ord_deriv,axis=1)*(1.0/dzdth[i,:,:])


    dudx = (dudxi*dydet-dudet*dydxi)/Jdet_xiet
    dvdy = (dvdet*dxdxi-dvdxi*dxdet)/Jdet_xiet
    dvdx = (dvdxi*dydet-dvdet*dydxi)/Jdet_xiet
    dudy = (dudet*dxdxi-dudxi*dxdet)/Jdet_xiet
    dwdx = (dwdxi*dydet-dwdet*dydxi)/Jdet_xiet
    dwdy = (dwdet*dxdxi-dwdxi*dxdet)/Jdet_xiet

    return  (dudx*dvdy+dudx*dwdz+dvdy*dwdz)-(dvdx*dudy+dudz*dwdx+dvdz*dwdy)
