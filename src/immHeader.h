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

#ifndef IMM_HEADER_H
#define IMM_HEADER_H

#include <stdint.h>

class IMMHeader
{

public:
    int32_t             mode;               //comp mode
    int32_t             compression;  //comp?
    unsigned char       date[32];        //date today
    unsigned char       prefix[16];       //?
    int32_t             number;     // filenum
    unsigned char       suffix[16];
    int32_t             monitor;    //0
    int32_t             shutter;    //0
    int32_t             row_beg;    //
    int32_t             row_end;      // whatever they are
    int32_t             col_beg;
    int32_t             col_end;
    int32_t             row_bin;      //1
    int32_t             col_bin;      //1
    int32_t             rows;         //
    int32_t             cols;
    int32_t             bytes;        //2
    int32_t             kinetics;     //0 part of ccd
    int32_t             kinwinsize;   //0
    double              elapsed;      //timestamp
    double              preset;       //exp time
    int32_t             topup;     //0
    int32_t             inject;       //0
    uint32_t            dlen;
    int32_t             roi_number;   //1

    int                 buffer_number;  //0
    unsigned int        systick;//0

    unsigned char       pv1[40];
    float               pv1VAL;
    unsigned char       pv2[40];
    float               pv2VAL;
    unsigned char       pv3[40];
    float               pv3VAL;
    unsigned char       pv4[40];
    float               pv4VAL;
    unsigned char       pv5[40];
    float               pv5VAL;
    unsigned char       pv6[40];
    float               pv6VAL;
    unsigned char       pv7[40];
    float               pv7VAL;
    unsigned char       pv8[40];
    float               pv8VAL;
    unsigned char       pv9[40];
    float               pv9VAL;
    unsigned char       pv10[40];
    float               pv10VAL;

    float               imageserver;  //make up
    float               CPUspeed;     //0

    enum {immver=12};

    int32_t             immversion;   //immver
    int                 corecotick;   //
    int                 cameratype;//160
    float               threshhold;   //my val

    //here is 632 bytes- or byte 0 to 631

    char byte632;//the is byte 632 counting from 0 to 632

    // 1024- (633+84 + 12)
    char empty_space[295];

    enum {z_len=84};
    enum {f_len = 12};

    //fill with 00's
    unsigned char ZZZZ[z_len];

    //fill with FF's
    unsigned char FFFF[f_len];

    enum {header_size=1024};
};

#endif