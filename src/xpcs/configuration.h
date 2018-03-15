/**

Copyright (c) 2016, UChicago Argonne, LLC. All rights reserved.

Copyright 2016. UChicago Argonne, LLC. This software was produced 
under U.S. Government contract DE-AC02-06CH11357 for Argonne National 
Laboratory (ANL), which is operated by UChicago Argonne, LLC for the 
U.S. Department of Energy. The U.S. Government has rights to use, 
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR 
UChicago Argonne, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR a
ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is 
modified to produce derivative works, such modified software should 
be clearly marked, so as not to confuse it with the version available 
from ANL.

Additionally, redistribution and use in source and binary forms, with 
or without modification, are permitted provided that the following 
conditions are met:

    * Redistributions of source code must retain the above copyright 
      notice, this list of conditions and the following disclaimer. 

    * Redistributions in binary form must reproduce the above copyright 
      notice, this list of conditions and the following disclaimer in 
      the documentation and/or other materials provided with the 
      distribution. 

    * Neither the name of UChicago Argonne, LLC, Argonne National 
      Laboratory, ANL, the U.S. Government, nor the names of its 
      contributors may be used to endorse or promote products derived 
      from this software without specific prior written permission. 

THIS SOFTWARE IS PROVIDED BY UChicago Argonne, LLC AND CONTRIBUTORS 
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL UChicago 
Argonne, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN 
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.

**/
#ifndef CONFIGURATION_
#define CONFIGURATION_

#include <string.h>
#include <map>
#include <vector>

#include "hdf5.h"

namespace xpcs  {

class Configuration  {

public:    
  static Configuration* instance();

  ~Configuration();

  int* getDQMap();
  int* getSQMap();
  int *PixelsPerStaticBin();

  int getFrameWidth();
  int getFrameHeight();

  int getFrameStartTodo();
  int getFrameEndTodo();
  int getFrameStart();
  int getFrameEnd();
  int getFrameTodoCount();
  int getFrameCount();
  int getDarkFrameStart();
  int getDarkFrameEnd();
  int getDarkFrames();

  void init(const std::string &path, const std::string &entry);

  std::map<int, std::map<int, std::vector<int>> > getBinMaps();

  int getTotalStaticPartitions();
  int getTotalDynamicPartitions();
  int getStaticWindowSize();

  std::string getFilename();
  std::string& getIMMFilePath();
  std::string& OutputPath();
  void setIMMFilePath(std::string& str);

  short* getPixelMask();
  int* getSbinMask();

  float getDetDpixX();
  float getDetDpixY();
  float getDetAdhuPhot();
  float getDetPreset();
  float getDetEfficiency();
  float getNormFactor();

  float getDarkThreshold();
  float getDarkSigma();

  bool getIsFlatFieldEnabled();
  double* getFlatField();

private:

  Configuration();

  std::string getString(const std::string &path);

  float getFloat(const std::string &path);

  int getInteger(const std::string &path);

  int* get2DTable(const std::string &path);

  double* get2DTableD(const std::string &path);

  void BuildQMap();
  
  double* flatfield;

  // Valid pixel mask - mark an entry as 1 in an array if the pixel is contained in any of the bins.
  short* m_validPixelMask;
  int* m_sbin;

  int m_totalStaticPartitions;
  int m_totalDynamicPartitions;

  // Map of dynamic bins to static bin to pixels.
  std::map<int, std::map<int, std::vector<int> >> m_mapping;

  int xdim;
  int ydim;
  int frameStart;
  int frameEnd;
  int frameStartTodo;
  int frameEndTodo;
  int darkFrameStart;
  int darkFrameEnd;
  int darkFrames;
  int m_staticWindow;

  int *pixels_per_bin;
  int *dqmap;
  int *sqmap;

  float m_detDpixX;
  float m_detDpixY;
  float m_detAdhupPhot;
  float m_detPreset;
  float m_detEfficiency;
  float m_normFactor;
  float darkThreshold;
  float darkSigma;

  hid_t file_id;

  // Flags for checkig if certain fields are enabled. 
  bool compression;
  bool flatfieldEnabled;

  std::string m_filename;
  std::string m_immFile;
  std::string output_path_;

  static Configuration *s_instance;
};

} // namespace xpcs

#endif