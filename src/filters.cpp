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

#include "filters.h"

void Filters::frameStride(SparseData *data, int strideSize)
{
    Configuration* conf = Configuration::instance();
    int w = conf->getFrameWidth();
    int h = conf->getFrameHeight();
    int frames = conf->getFrameTodoCount();
    int pixels = w * h;

    std::vector<int> validPixels = data->getValidPixels();

    // #pragma omp parallel for default(none) shared(validPixels, delays_per_level, frames, pixels, G2s, IFs, IPs, data)

    for (int i = 0; i < validPixels.size(); i++)
    {
        Row *row = data->get(validPixels.at(i));
        std::vector<int> indxPtr = row->indxPtr;
        std::vector<float> valPtr = row->valPtr;

        int fstart = 0;
        int fend = fstart + strideSize;
        int findex = 0;

        for (int x = 0 ; x < indxPtr.size() ; x++)
        {
            int ind = indxPtr[x];
            float val = valPtr[x];

            if ( (ind % strideSize) == 0)
            {
                indxPtr[findex] = ind;
            }

        }

    }

}

void Filters::frameAverage(SparseData *data, int avgSize)
{
    Configuration* conf = Configuration::instance();
    int w = conf->getFrameWidth();
    int h = conf->getFrameHeight();
    int frames = conf->getFrameTodoCount();
    int pixels = w * h;

    std::vector<int> validPixels = data->getValidPixels();

    // #pragma omp parallel for default(none) shared(validPixels, delays_per_level, frames, pixels, G2s, IFs, IPs, data)
    for (int i = 0; i < validPixels.size(); i++)
    {
        Row *row = data->get(validPixels.at(i));
        std::vector<int> indxPtr = row->indxPtr;
        std::vector<float> valPtr = row->valPtr;

        int fstart = 0;
        int fend = fstart + strideSize;
        int findex = 0;

        for (int x = 0 ; x < indxPtr.size() ; x++)
        {
            int ind = indxPtr[x];
            float val = valPtr[x];

            if ( (ind % strideSize) == 0)
            {
                indxPtr[findex] = ind;
            }

        }

    }
}