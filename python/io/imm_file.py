import struct
import numpy as np

class Header(dict):
  def __init__(self, fp):
    self._dict = {}
    self.isValid = True
    self.__helper(fp)

  def update(self, fp):
    self.__helper(fp)

  def __helper(self, f):
    
    # TODO: Replace it with a better EOF check
    if f.read(4) == '':
      self.isValid = False
      return

    compression = struct.unpack('i', f.read(4))[0]
    f.seek(108, 1) # 116
    size = struct.unpack('i', f.read(4))[0]
    f.seek(32, 1) #152
    dlen = struct.unpack('i', f.read(4))[0]
    f.seek(4, 1) #160
    bufferno = struct.unpack('i', f.read(4))[0]

    self['compression'] = compression == 6
    self['bytes_per_pixel'] = size
    self['pixels'] = dlen
    self['frame_no'] = bufferno

    # Move it to the beginning of the data. 
    f.seek((1024 - 164), 1)

    #TODO Need additional attributes. 

class IMM:
  def __init__(self, filename, mode="rb"):
    self.imm_file = filename
    self.compression = None
    self.filePtr = None # File byte pointer
    self.frameIndex = 0

    # TODO:
    # Check mode string against supported values.

    try:
      self.filePtr = open(self.imm_file, mode)
    except Exception as e:
      print("Failed to read imm file %s"%self.imm_file)
      self.filePtr = None

    try:
      self.header = Header(self.filePtr)
    except Exception as ee:
      print ee
      print("IMM file doesn't seems to be of right type")

  def array(self, startIndex=0, endIndex=-1):
    
    result_array = np.empty((266, self.header["pixels"]), np.uint16)

    while self.header.isValid:
      compression = self.header["compression"]
      pixels = self.header["pixels"]
      bytes_per_pixel = self.header["bytes_per_pixel"]
      frameIndex = self.header["frame_no"]

      #print self.filePtr.tell()
      r = np.fromfile(self.filePtr, dtype=np.uint16, count=pixels)
      result_array[self.header["frame_no"]] = r

      #self.filePtr.seek(pixels * bytes_per_pixel, 1)
      self.header.update(self.filePtr)

    return result_array
