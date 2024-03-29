import numpy as np
import os
import glob
#
# define some custom parameters, not defined in the DNS code
#
iseek      = 0              # number of bytes to skip relative to the origin of the binary file (0 for SNaC)
iprecision = 8              # precision of real-valued data
if(    iprecision == 4):
    my_dtype = 'float32'
else:
    my_dtype = 'float64'
r0_g = np.array([0.,0.,0.]) # domain origin
non_uniform_grid = True
iskip      = 1
#
# retrieve number of blocks from the number of geo files
#
geofiles  = "geometry_b_???.out"
nblocks = np.size(glob.glob(geofiles))
class block:
    def __init__(self,iblock,blockname,saves,nsaves,nflds,gridfiles,n):
        self.id   = iblock
        self.name = blockname
        self.saves = saves
        self.nsaves    = nsaves
        self.nflds     = nflds
        self.gridfiles = gridfiles
        self.n         = n
blocks = []
block_suffix = '_b_'
logfile_prefix = input("Name of the log file written by SNaC (excluding the block-specific suffix) [log_visu_3d]: ") or "log_visu_3d"
grid_prefix = input("Name to be appended to the grid files to prevent overwriting [""]: ") or ""
for iblock in range(1,nblocks+1):
    blockname = block_suffix+"{:3}".format(str(iblock).zfill(3))
    #
    # define data type and
    # read saved data log
    #
    dtype_saves = np.dtype([                                                \
                            ('file' , 'U100'), ('variable', 'U100'),        \
                            ('imin' , int), ('jmin' , int), ('kmin' , int), \
                            ('imax' , int), ('jmax' , int), ('kmax' , int), \
                            ('istep', int), ('jstep', int), ('kstep', int), \
                            ('time', float), ('isave', int)                 \
                           ])
    geofile  = "geometry"+blockname+".out"
    logfile  = logfile_prefix+blockname+".out"
    gridname = grid_prefix+blockname
    xgridfile = "x"+gridname+'.bin'
    ygridfile = "y"+gridname+'.bin'
    zgridfile = "z"+gridname+'.bin'
    saves = np.loadtxt(logfile, dtype=dtype_saves)
    #
    # remove duplicates
    #
    saves = np.unique(saves)
    #
    # sort elements by increasing isave
    #
    saves = np.sort(saves, order='isave')
    #
    # harvest some information from the log file
    #
    nelements = saves.size
    nflds     = 0
    isave = saves['isave'][0]
    while(isave == saves['isave'][nflds] and nflds < nelements-1): nflds += 1
    if(nflds == nelements-1): nflds += 1
    nsaves = int(nelements/nflds)
    nmin  = np.array([saves['imin' ][0], saves['jmin' ][0], saves['kmin' ][0]])
    nmax  = np.array([saves['imax' ][0], saves['jmax' ][0], saves['kmax' ][0]])
    nstep = np.array([saves['istep'][0], saves['jstep'][0], saves['kstep'][0]])
    n = ((nmax-nmin+1)/nstep).astype(int)
    #
    # retrieve some computational parameters
    #
    data = np.loadtxt(geofile, comments = "!", max_rows = 4)
    ng = (data[1,:]-data[0,:]+1).astype('int')
    l  = data[3,:]-data[2,:]
    dl = l/(1.*ng)
    r0 = data[2,:]
    #
    # generate grid files
    #
    x = np.arange(r0[0]+dl[0]/2.,r0[0]+l[0],dl[0])
    y = np.arange(r0[1]+dl[1]/2.,r0[1]+l[1],dl[1])
    z = np.arange(r0[2]+dl[2]/2.,r0[2]+l[2],dl[2])
    if os.path.exists(xgridfile): os.remove(xgridfile)
    if os.path.exists(ygridfile): os.remove(ygridfile)
    if os.path.exists(zgridfile): os.remove(zgridfile)
    if(non_uniform_grid):
        f   = open('grid_x'+blockname+'.bin','rb')
        grid_x = np.fromfile(f,dtype=my_dtype)
        f.close()
        f   = open('grid_y'+blockname+'.bin','rb')
        grid_y = np.fromfile(f,dtype=my_dtype)
        f.close()
        f   = open('grid_z'+blockname+'.bin','rb')
        grid_z = np.fromfile(f,dtype=my_dtype)
        f.close()
        grid_x = np.reshape(grid_x,(ng[0],4),order='F')
        grid_y = np.reshape(grid_y,(ng[1],4),order='F')
        grid_z = np.reshape(grid_z,(ng[2],4),order='F')
        x = r0_g[0] + grid_x[:,1]
        y = r0_g[1] + grid_y[:,1]
        z = r0_g[2] + grid_z[:,1]
    x[nmin[0]-1:nmax[0]:nstep[0]].astype(my_dtype).tofile(xgridfile)
    y[nmin[1]-1:nmax[1]:nstep[1]].astype(my_dtype).tofile(ygridfile)
    z[nmin[2]-1:nmax[2]:nstep[2]].astype(my_dtype).tofile(zgridfile)
    blocks.append(block(iblock,blockname,saves,nsaves,nflds,np.array([xgridfile,ygridfile,zgridfile]),n))
    print('Block # {:3} log files parsed and grid files generated.'.format(str(iblock).zfill(3)))
#
# small sanity check
#
same_temporal_series = True
for iblock in range(nblocks):
    same_temporal_series = same_temporal_series and \
                           blocks[iblock].nsaves == blocks[iblock-1].nsaves and \
                           blocks[iblock].nflds == blocks[iblock-1].nflds and \
                           all(blocks[iblock].saves["time"][:]==blocks[iblock-1].saves["time"][:])
import sys
if(not same_temporal_series): sys.exit("Error: all block files must contain the same time series")
nsaves = blocks[0].nsaves
times  = blocks[0].saves["time"][:]
nflds  = blocks[0].nflds
#
# write xml file
#
from xml.etree import ElementTree
from xml.dom import minidom
from xml.etree.ElementTree import Element, SubElement, Comment
Xdmf = Element("Xdmf", attrib = {"xmlns:xi": "http://www.w3.org/2001/XInclude", "Version": "2.0"})
domain = SubElement(Xdmf, "Domain")
for iblock in range(nblocks):
    blockname = blocks[iblock].name
    n = blocks[iblock].n
    topology = SubElement(domain, "Topology", attrib = {"name": "TOPO"+blockname, "TopologyType": "3DRectMesh", "Dimensions" : "{} {} {}".format(n[2], n[1], n[0])})
    geometry = SubElement(domain, "Geometry", attrib = {"name": "GEO" +blockname, "GeometryType": "VXVYVZ"})
    dataitem = SubElement(geometry, "DataItem", attrib = {"Format": "Binary", "DataType": "Float", "Precision": "{}".format(iprecision), "Endian": "Native", "Dimensions": "{}".format(n[0])})
    dataitem.text = blocks[iblock].gridfiles[0]
    dataitem = SubElement(geometry, "DataItem", attrib = {"Format": "Binary", "DataType": "Float", "Precision": "{}".format(iprecision), "Endian": "Native", "Dimensions": "{}".format(n[1])})
    dataitem.text = blocks[iblock].gridfiles[1]
    dataitem = SubElement(geometry, "DataItem", attrib = {"Format": "Binary", "DataType": "Float", "Precision": "{}".format(iprecision), "Endian": "Native", "Dimensions": "{}".format(n[2])})
    dataitem.text = blocks[iblock].gridfiles[2]
grid = SubElement(domain, "Grid", attrib = {"Name": "TimeSeries", "GridType": "Collection",  "CollectionType": "Temporal"})
#time = SubElement(grid, "Time", attrib = {"TimeType":"List"})
#dataitem = SubElement(time, "DataItem", attrib = {"Format": "XML", "NumberType": "Float", "Dimensions": "{}".format(nsaves)})
#dataitem.text = ""
#for ii in range(nsaves):
    #dataitem.text += "{:15.6E}".format(times[ii*nflds]) + " "
istart = 0
iend   = nsaves
iskip  = iskip
for ii in range(istart, iend, iskip):
    grid_blk = SubElement(grid, "Grid", attrib = {"Name": "T{:7}".format(str(blocks[iblock-1].saves['isave'][ii*nflds]).zfill(7))+"_ALL", "GridType": "Collection", "CollectionType": "Spatial"})
    time_blk = SubElement(grid_blk, "Time", attrib = {"Value": "{:15.6E}".format(times[ii*nflds])})
    for iblock in range(nblocks):
        blockname = blocks[iblock].name
        n = blocks[iblock].n
        grid_fld = SubElement(grid_blk, "Grid", attrib = {"Name": "T{:7}".format(str(blocks[iblock-1].saves['isave'][ii*nflds]).zfill(7))+blockname, "GridType": "Uniform"})
        topology = SubElement(grid_fld, "Topology", attrib = {"Reference": "/Xdmf/Domain/Topology[{}]".format(iblock+1)})
        geometry = SubElement(grid_fld, "Geometry", attrib = {"Reference": "/Xdmf/Domain/Geometry[{}]".format(iblock+1)})
        for jj in range(nflds):
            index = ii*nflds+jj
            #
            # if vector, skip second and third components from the loop, and write the three files at once
            #
            if( (blocks[iblock].saves['variable'][index-1+0].endswith('_X') and blocks[iblock].saves['variable'][index-1+1].endswith('_Y') and blocks[iblock].saves['variable'][index-1+2].endswith('_Z') and \
                 blocks[iblock].saves['variable'][index-1+0][0:-2] ==           blocks[iblock].saves['variable'][index-1+1][0:-2] ==           blocks[iblock].saves['variable'][index-1+2][0:-2]) or \
                (blocks[iblock].saves['variable'][index-2+0].endswith('_X') and blocks[iblock].saves['variable'][index-2+1].endswith('_Y') and blocks[iblock].saves['variable'][index-2+2].endswith('_Z') and \
                 blocks[iblock].saves['variable'][index-2+0][0:-2] ==           blocks[iblock].saves['variable'][index-2+1][0:-2] ==           blocks[iblock].saves['variable'][index-2+2][0:-2]) \
              ): continue
            #
            # vector
            #
            if(blocks[iblock].saves['variable'][index+0].endswith('_X') and blocks[iblock].saves['variable'][index+1].endswith('_Y') and blocks[iblock].saves['variable'][index+2].endswith('_Z') and \
               blocks[iblock].saves['variable'][index+0][0:-2] ==           blocks[iblock].saves['variable'][index+1][0:-2] ==           blocks[iblock].saves['variable'][index+2][0:-2]):
               attribute = SubElement(grid_fld, "Attribute", attrib = {"AttributeType": "Vector", "Name": "{}".format(blocks[iblock].saves['variable'][index][:-2]), "Center": "Node"})
               attribute = SubElement(attribute, "DataItem", attrib = {"Function": "JOIN($0, $1, $2)", "ItemType": "Function", "Dimensions": "{} {} {} {}".format(n[2], n[1], n[0], 3)})
               for q in range(3):
                  dataitem = SubElement(attribute, "DataItem", attrib = {"Format": "Binary", "DataType": "Float", "Precision": "{}".format(iprecision), "Endian": "Native", "Seek": "{}".format(iseek), "Dimensions": "{} {} {}".format(n[2], n[1], n[0])})
                  dataitem.text = blocks[iblock].saves['file'][index+q]
            #
            # scalar
            #
            else:
               attribute = SubElement(grid_fld, "Attribute", attrib = {"AttributeType": "Scalar", "Name": "{}".format(blocks[iblock].saves['variable'][index]), "Center": "Node"})
               dataitem = SubElement(attribute, "DataItem", attrib = {"Format": "Binary", "DataType": "Float", "Precision": "{}".format(iprecision), "Endian": "Native", "Seek": "{}".format(iseek), "Dimensions": "{} {} {}".format(n[2], n[1], n[0])})
               dataitem.text = blocks[iblock].saves['file'][index]
output = ElementTree.tostring(Xdmf, 'utf-8')
output = minidom.parseString(output)
output = output.toprettyxml(indent="    ",newl='\n')
#
# write visualization file
#
outfile = input("Name of the output file [viewfld_DNS.xmf]: ") or "viewfld_DNS.xmf"
xdmf_file = open(outfile, 'w')
xdmf_file.write(output)
xdmf_file.close
#
# workaround to add the DOCTYPE line
#
xdmf_file = open(outfile, "r")
contents = xdmf_file.readlines()
xdmf_file.close()
contents.insert(1, '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
xdmf_file = open(outfile, "w")
contents = "".join(contents)
xdmf_file.write(contents)
xdmf_file.close()
