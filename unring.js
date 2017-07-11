/*
Copyright 2017 Paul Reggentin

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


*/

math = require('./kmath.js');
//require("./complex.js");
fs = require("fs");
require('./nifti.js');
var FFT = require("fft");



//RUN FROM COMMAND LINE WITH NODE
/*
    var path_in = process.argv[2];
    var path_out = process.argv[3];


run_unring(path_in,path_out,30,1,4,function(){});

console.log('done!!');
*/



//RUN AS SCRIPT


//set filenames and process input nifti file
var filename = 't2zoom_test.nii'

var outfile = 't2zoom_unring.nii'

//var raw = fs.readFileSync(filename);

//var nii = parse(raw);

//perform the unringing
run_unring(filename,outfile,30,1,4,function(){});



function run_unring(filename,outfile,nsh,minW,maxW,done)
///(img_in, is_cplx, img_out, dim, nsh, minW,maxW)
{
    var raw = fs.readFileSync(filename);
    var nii = parse(raw.buffer);
    
    var data = nii.data;
    var sz = nii.sizes;

    exportnii(nii.data,'before.nii');

    // happy unringing
    var raw_slice = new Float64Array(sz[0]*sz[1]);
    var unrung_slice = new Float64Array(sz[0]*sz[1]*2);
    var unrung_int16 = new Int16Array(sz[0]*sz[1]);
    
    for(var slice = 0; slice < sz[2]; slice++){
        var st = performance.now();
        var offset=slice*nii.sizes[0]*nii.sizes[1];
        
        raw_slice = nii.data.slice(offset,offset+nii.sizes[0]*nii.sizes[1]);
        
        unring2d(raw_slice,false,unrung_slice,sz,nsh,minW,maxW);
        
        //convert unrung data to int16 for output
        for(var i = 0; i < unrung_int16.length; i++){
            nii.data[i+offset] = unrung_slice[2*i] + .5;
        }
        var en = performance.now();
        console.log('total time:',en-st);
        console.log(slice);
        
    }//for slice
    exportnii(nii.data,'after.nii');
    fs.writeFileSync(outfile, raw);
    done();
}



/*
    img_in: 2D image, with possible ringing, stored as 1D Float64Array 
    img_out: 2D image stored as 1D Float64Array which will store the unrung image
    dim: 2-elt array with the dimensions of the input image
    nsh: number of shifts to try
    minW: minimum of the shift window
    maxW: maximum of the shift window
*/
function unring2d(img_in, is_cplx, img_out, dim, nsh, minW,maxW){

    //create the fft variables for dim 1 and 2 of the image
    fft1 = FFT.complex(dim[0],false);
    ifft1 = FFT.complex(dim[0],true);
    fft2 = FFT.complex(dim[1],false);
    ifft2 = FFT.complex(dim[1],true);

    //create the fourier transforms of the image and its transpose
    var fac = 2;
    if(is_cplx){
        fac = 1;
    }
    var img_in_ft = new Float64Array(img_in.length*fac);
    var img_in_t_ft = new Float64Array(img_in.length*fac);

    fft_2d(img_in,dim,img_in_ft,is_cplx);
    transpose_cplx(img_in_ft,dim,img_in_t_ft);
    //(the 2d fft of the transpose of X i the hermitian transpose of the 2d fft of x)

    //create and apply saddle filter
    for(var j = 0; j < dim[1]; j++){
        var cj = (1+Math.cos(2*Math.PI*j/dim[1]))*0.5;
        for(var i = 0; i < dim[0]; i++){
            var ci = (1+Math.cos(2*Math.PI*i/dim[0]))*0.5;           
            var eps = .000000000000001
            var scale = 1/(ci+cj+eps);//the denominator of the filter value, eps to avoid /0 error
            img_in_ft[2*(i+j*dim[0])] = img_in_ft[2*(i+j*dim[0])] * cj * scale;
            img_in_ft[2*(i+j*dim[0]) + 1] = img_in_ft[2*(i+j*dim[0]) + 1] * cj * scale;
            img_in_t_ft[2*(j + i*dim[1])] = img_in_t_ft[2*(j + i*dim[1])] * ci * scale;  
            img_in_t_ft[2*(j + i*dim[1]) + 1] = img_in_t_ft[2*(j + i*dim[1]) + 1] * ci * scale;  
        }//for i
    }//for j

    //convert the filtered images back to spatial domain
    var filt1 = new Float64Array(img_in.length*2);
    var filt2 = new Float64Array(img_in.length*2);
    ifft_2d(img_in_ft,dim,filt1);
    ifft_2d(img_in_t_ft,dim_t,filt2);

    //unring the filtered images
    unring1d(filt1,dim[0],dim[1],nsh,minW,maxW,function(){});
    unring1d(filt2,dim[1],dim[0],nsh,minW,maxW,function(){});

    //add the two filtered images together and store them in img_out
    for(var i = 0; i < dim[0]; i++){
        for(var j = 0; j < dim[1]; j++){
            img_out[2*(i+j*dim[0])] = filt1[2*(i+j*dim[0])] + filt2[2*(j + i*dim[1])];
            img_out[2*(i+j*dim[0]) + 1] = filt1[2*(i+j*dim[0]) + 1] + filt2[2*(j + i*dim[1]) + 1];
        }
    }
}//unring2d

/*
    transpose a real-valued 2-d image 
    img: 1-D vector storing the original image
    dim: the dimensions of the 2-d image stored in img vector
    img_t: 1-D vector that will hold the transposed image 
*/
function transpose_re(img,dim,img_t){
    for(var i = 0; i < dim[0]; i++){
        for(var j = 0; j < dim[1]; j++){
            img_t[j+i*dim[1]] = img[i+j*dim[0]];
        }
    }
}

/*
    transpose a real-valued 2-d image 
    img: 1-D vector storing the original image 
            in the form [real0 imag0 real1 imag1 ... ]
    dim: the dimensions of the 2-d image stored in img vector
    img_t: 1-D vector that will hold the transposed image 
*/
function transpose_cplx(img,dim,img_t){
    for(var i = 0; i < dim[0]; i++){
        for(var j = 0; j < dim[1]; j++){

            img_t[2*(j+i*dim[1])] = img[2*(i+j*dim[0])];
            img_t[2*(j+i*dim[1]) + 1] = img[2*(i+j*dim[0]) + 1];
            

        }
    }
}

/*
    insert a source vector into a destination vector
    src: the vector to be inserted
    dest: the vector that src is inserted into
    n: number of elements of src to insert
    offset: the position in dest that src is inserted into
*/
function copy_into(src,dest,n,offset){
    for(var i = 0; i < n; i++){
        dest[offset+i] = src[i];
    }
}

/*
    perform a 2-dimensional fourier transform 
    img: the image to be transformed, stored as a 1D vector. 
        can be real or complex, stored as [real0 imag0 real1 imag1 ...]
    dim: the dimensions of the 2D image that 'img' represents
    img_ft: the vector that will hold the fourier transform of img
    is_cplx: boolean flag for type of img
*/
function fft_2d(img, dim, img_ft, is_cplx){
    //create fft objects
    fft0 = new FFT.complex(dim[0],false);
    fft1 = new FFT.complex(dim[1],false);

    //create vectors that will store intermediate steps
    var frow = new Float64Array(2*dim[0]);
    var intermed = new Float64Array(2*dim[0]*dim[1]);
    var intermed_t = new Float64Array(2*dim[0]*dim[1]);
    var fft_trans = new Float64Array(2*dim[0]*dim[1]);

    //assign the correct input row length
    var fft_type;
    var row_len;
    if(is_cplx){
        fft_type = 'complex';
        row_len = 2*dim[0];
    } else {
        fft_type = 'real';
        row_len = dim[0];
    }

    //perform 1D fft on each of the rows of img
    for(var i = 0; i < dim[1]; i++){
        fft0.simple(frow,img.slice(i*row_len,(i+1)*row_len),fft_type);
        copy_into(frow,intermed,2*dim[0],i*2*dim[0]);
    }

    //transpose the row transformed array
    transpose_cplx(intermed,dim,intermed_t);

    //fft each of the rows (which are columns from the original image)
    var col_len = 2*dim[1];
    var fcol = new Float64Array(col_len);
    for(var i = 0; i < dim[0]; i++){
        fft1.simple(fcol,intermed_t.slice(i*col_len,(i+1)*col_len),'complex');
        copy_into(fcol,fft_trans,col_len,i*col_len);
    }

    //transpose back to original dimentions
    dim_t = new Float64Array([dim[1], dim[0]]);
    transpose_cplx(fft_trans,dim_t,img_ft);
}

//inverse fourier transform
//same as fft_2d but all of the 1D ffts are inverse DFTs
function ifft_2d(kspace, dim, img_space){
    fft0 = new FFT.complex(dim[0],true);
    fft1 = new FFT.complex(dim[1],true);

    //fft each of the rows
    var frow = new Float64Array(2*dim[0]);
    var intermed = new Float64Array(2*dim[0]*dim[1]);
    var intermed_t = new Float64Array(2*dim[0]*dim[1]);
    var fft_trans = new Float64Array(2*dim[0]*dim[1]);

    var fft_type = 'complex';
    var row_len = 2*dim[0];

    for(var i = 0; i < dim[1]; i++){

        fft0.simple(frow,kspace.slice(i*row_len,(i+1)*row_len),fft_type);
        copy_into(frow,intermed,2*dim[0],i*2*dim[0]);
        
    }

    //transpose the row transformed
    transpose_cplx(intermed,dim,intermed_t);

    //fft each of the rows (which are columns from the original image)
    var col_len = 2*dim[1];
    var fcol = new Float64Array(col_len);
    for(var i = 0; i < dim[0]; i++){
        fft1.simple(fcol,intermed_t.slice(i*col_len,(i+1)*col_len),'complex');
        copy_into(fcol,fft_trans,col_len,i*col_len);
    }
    //transpose back to original dimentions
    dim_t = new Float64Array([dim[1], dim[0]]);
    transpose_cplx(fft_trans,dim_t,img_space);
    
    //perform scaling
    scale = 1/(dim[0]*dim[1]);
    for(var i = 0; i < img_space.length; i++){
        img_space[i] = img_space[i]*scale; 
    }
}



/*
perform the unringing algorithm on each row of a 2-d array
ARGS: 
    data     - the numlines-by-n image array
    n        - number of pixels per row
    numlines - number of rows to analzye
    nsh      - number of subvoxel shift values to try in one direction (total # shifts will be 2*nsh+1)
    minW     - min of window
    maxW     - max of window
    callback - callback function
*/
function unring1d(data, n, numlines, nsh, minW, maxW, callback){

    //DEBUG: time vars
    /*
    var end_time;
    var start_time;
    var shift_time= 0;
    var total_time= 0;
    var calc_time= 0;
    */
    
   
    //create the length variables
    var totalShifts = 2*nsh+1;
    var totalSamps = 2*n;
    
    //create the necessary variables
    var fft = new FFT.complex(n,false);
    var ifft = new FFT.complex(n,true);

    //this variable holds the frequencydomain version of a row, so that it can be shifted
    var rowFreq = new Float64Array(totalSamps);
    var rowFreqShift = new Float64Array(totalSamps);
    var rowTimeShift = new Float64Array(totalSamps);

    //these variables hold the rows shifted in time and freq domains
    var timeShifts = new Float64Array(totalSamps*totalShifts);
    //var freqShifts = new Float64Array(totalSamps*totalShifts);

    //these arrays hold the measure of ringing at successive shifts to the right and leftDiffs
    //ringing is measured as the sum of the pairwise differences of adjacent pixels
    var rightDiffs = new Float64Array(totalShifts);
    var leftDiffs = new Float64Array(totalShifts);


    //create the shift array
    var shifts = new Float64Array(totalShifts);
    shifts[0] = 0;
    for(var i = 0; i < nsh; i++){
        shifts[i+1] = (i+1);
        shifts[i+1+nsh] = 0-(i+1);
    }
    
    //the angle corresponding to a certain linear shift
    var phi = new Float64Array(2);
    //the angle that increases for each successive frequencz
    var ang = new Float64Array(2);


    //var maxn;
    if(n%2 == 1) {
        var maxn = (n-1)/2;
    } else {
        var maxn = n/2-1;
    }
    
    //generate the frequency ramp
    var freqRamps = new Float64Array(2*maxn*shifts.length);
    var ramplen = 2*maxn;
    for(s = 1; s < shifts.length; s++){
        var p =  shifts[s] * Math.PI / (n*nsh);  
        phi[0] = Math.cos(p);
        phi[1] = Math.sin(p);
        ang[0] = 1;
        ang[1] = 0;
        for(var i =0; i < maxn; i++){
            var tmp = ang[0];
            ang[0] = phi[0]*ang[0] - phi[1]*ang[1];
            ang[1] = tmp*phi[1] + phi[0]*ang[1];
            freqRamps[s*ramplen + 2*i] = ang[0];
            freqRamps[s*ramplen + 2*i + 1] = ang[1];
        }
    }

    //apply shift in frequency domain, then in time domain
    for (var line = 0; line < numlines; line++){
        var line_idx = line*totalSamps;
         //DEBUG
        //start_time = performance.now();
        for(var i = 0; i < totalSamps; i++){
            timeShifts[i] = data[line_idx+i];
        }
        
        

        
        fft.simple(rowFreq,data.slice(line_idx,line_idx+totalSamps),'complex');
        
        //shift the frequency data by each value in shifts
        //s is shift number
        //var angle_scale = Math.PI /(n*nsh);
        for(var s = 1; s < shifts.length; s++){


            /*
            //set the angle values
            var p =  shifts[s] * Math.PI /(n*nsh);
            phi[0] = Math.cos(p);
            phi[1] = Math.sin(p);
            ang[0] = 1;
            ang[1] = 0;
            */
            
            //set the dc term and nyquist frequency
            rowFreqShift[0] = rowFreq[0];
            rowFreqShift[1] = rowFreq[1];
            
            //if the Nquist frequency is included, set it to 0
            if(n%2==0) {
                rowFreqShift[n] = 0;
                rowFreqShift[ n + 1] = 0;
            }
            
            var tmp;
            for(var i =0; i < maxn; i++){
                var L = (i+1);
                var FR = 2*i;
                /*
                var tmp = ang[0];
                ang[0] = phi[0]*ang[0] - phi[1]*ang[1];
                ang[1] = tmp*phi[1] + phi[0]*ang[1];
                if(Math.abs(ang[0]-freqRamps[s*ramplen + FR]) > 0.00000001 || Math.abs(ang[1]-freqRamps[s*ramplen + FR + 1])>0.00000001){
                    console.log('wrong!');
                }
                */


                rowFreqShift[2*L] = freqRamps[s*ramplen + FR]*rowFreq[2*L] - freqRamps[s*ramplen + FR + 1]*rowFreq[2*L+1];
                rowFreqShift[2*L+1] = freqRamps[s*ramplen + FR]*rowFreq[2*L+1] + freqRamps[s*ramplen + FR + 1]*rowFreq[2*L];
                
                //rowFreqShift[2*L] = ang[0]*rowFreq[2*L] - ang[1]*rowFreq[2*L+1];
                //rowFreqShift[2*L+1] = ang[0]*rowFreq[2*L+1] + ang[1]*rowFreq[2*L];
                


                L = (n-1-i);
                rowFreqShift[2*L] = freqRamps[s*ramplen + FR]*rowFreq[2*L] + freqRamps[s*ramplen + FR + 1]*rowFreq[2*L+1];
                rowFreqShift[2*L+1] = freqRamps[s*ramplen + FR]*rowFreq[2*L+1] - freqRamps[s*ramplen + FR + 1]*rowFreq[2*L];
                
                
                //rowFreqShift[2*L] = ang[0]*rowFreq[2*L] + ang[1]*rowFreq[2*L+1];
                //rowFreqShift[2*L+1] = ang[0]*rowFreq[2*L+1] - ang[1]*rowFreq[2*L];
                
            }//for i (freq bin to shift)
            


            //perform the inverse fft to get signal shifted in time domain
            ifft.simple(rowTimeShift,rowFreqShift,'complex');

            
            //copy the time shifted array into the dictionary of shifted values
            //move this to the loop above
            for(var i = 0; i < totalSamps; i++){
                timeShifts[s*totalSamps + i] = rowTimeShift[i] / n;
            }
        }//for s (shift number)

        


        //sum the adjacent differences to get the initial ringing metrics for pixel 0
        for(var s = 0; s<totalShifts; ++s){
            var offset = s*totalSamps;
            rightDiffs[s] = 0;
            leftDiffs[s] = 0;
            //var pix = 0;
            //sum the differences in the window for 
            for(var d = minW; d <= maxW; d++){
                
                var rightIdx1 = (2*(d) + totalSamps)%totalSamps;
                var rightIdx2 = (rightIdx1+2)%totalSamps;
                var t = -2;
                var leftIdx1 = (t*(d) + totalSamps)%totalSamps;
                var leftIdx2 = (leftIdx1-2)%totalSamps;
                
                //manhattan distance
                
                rightDiffs[s] += Math.abs( timeShifts[offset + rightIdx1] - timeShifts[offset + rightIdx2]);
                //rightDiffs[s] += Math.abs( timeShifts[offset + rightIdx1 + 1] - timeShifts[offset + rightIdx2 + 1]);
            
                leftDiffs[s] += Math.abs( timeShifts[offset + leftIdx1] - timeShifts[offset + leftIdx2]);
                //leftDiffs[s] += Math.abs( timeShifts[offset + leftIdx1 + 1] - timeShifts[offset + leftIdx2 + 1]);
                
                var diff = timeShifts[offset + rightIdx1] - timeShifts[offset + rightIdx2];
                if (diff>0) {
                    rightDiffs[s] += diff;
                } else {
                    rightDiffs[s] -=diff;
                }

                diff = timeShifts[offset + leftIdx1] - timeShifts[offset + leftIdx2];
                if (diff>0) {
                    leftDiffs[s] += diff;
                } else {
                    leftDiffs[s] -= diff;
                }

            }//for d

        }//for s
        
        //DEBUG
        //end_time = performance.now();
       //shift_time += end_time-start_time;
        


        //for each pixel find the minimal ringing measure & shift to that ringing measure
        //DEBUG
         //start_time = performance.now();
         
        for(var pix = 0; pix < n; pix++){
            var minDiff = 99999999999;
            var minIndex = 0;


            
        
         for(var s = 0; s< totalShifts; s++){

                var offset = s*totalSamps;
                
                //console.log(rightDiffs);
                
                if (rightDiffs[s] < minDiff){
                    minDiff = rightDiffs[s];
                    minIndex = s;
                } if (leftDiffs[s] < minDiff){
                    minDiff = leftDiffs[s];
                    minIndex = s;
                }

                //get the ringing measure for the successive pixel by removing the distance measure on one end and adding the distance measure on the other end.
                var oldRight1 = (2*(pix+minW) + totalSamps) % totalSamps;
                var oldRight2 = (oldRight1 + 2) % totalSamps;
                var newRight1 = (2*(pix+maxW+1) + totalSamps) % totalSamps;
                var newRight2 = (newRight1 + 2) % totalSamps;

                var oldLeft1 = (2*(pix - maxW) + totalSamps)%totalSamps;
                var oldLeft2 = (oldLeft1-2)%totalSamps;
                var newLeft1 = (2*(pix - minW - 1) + totalSamps)%totalSamps;
                var newLeft2 = (newLeft1-2)%totalSamps;
                
                //manhattan distance
                //subtract off old bounds
                
                //rightDiffs[s] -= Math.abs( timeShifts[offset + oldRight1] - timeShifts[offset + oldRight2]);
                //rightDiffs[s] -= Math.abs( timeShifts[offset + oldRight1 + 1] - timeShifts[offset + oldRight2 + 1]);
                //leftDiffs[s] -= Math.abs( timeShifts[offset + oldLeft1] - timeShifts[offset + oldLeft2]);
                //leftDiffs[s] -= Math.abs( timeShifts[offset + oldLeft1 + 1] - timeShifts[offset + oldLeft2 + 1]);
                
                //add on new bounds
                //rightDiffs[s] += Math.abs( timeShifts[offset + newRight1] - timeShifts[offset + newRight2]);
                //rightDiffs[s] += Math.abs( timeShifts[offset + newRight1 + 1] - timeShifts[offset + newRight2 + 1]);
                //leftDiffs[s]  += Math.abs( timeShifts[offset + newLeft1]  - timeShifts[offset + newLeft2]);
                //leftDiffs[s] += Math.abs( timeShifts[offset + newLeft1 + 1] - timeShifts[offset + newLeft2 + 1]);
                ////console.log(Math.abs( timeShifts[offset + newLeft1 + 1] - timeShifts[offset + newLeft2 + 1]));
                
                
                //custom abs function
                var diff = timeShifts[offset + oldRight1] - timeShifts[offset + oldRight2];
                if(diff > 0) {
                    rightDiffs[s] -= diff;
                } else {
                    rightDiffs[s] += diff;
                }

                diff = timeShifts[offset + oldLeft1] - timeShifts[offset + oldLeft2];
                if(diff > 0) {
                    leftDiffs[s] -= diff;
                } else {
                    leftDiffs[s] += diff;
                }

                diff = timeShifts[offset + newRight1] - timeShifts[offset + newRight2];
                if(diff > 0){
                    rightDiffs[s] += diff;
                } else {
                    rightDiffs[s] -= diff;
                }

                diff = timeShifts[offset + newLeft1]  - timeShifts[offset + newLeft2];
                if(diff > 0){
                    leftDiffs[s] += diff;
                } else {
                    leftDiffs[s] -= diff;
                }
            }//for s
            
            


            
            //perform ye olde reinterpolation
            /*
            var dif0re = timeShifts[minIndex*totalSamps + (2*(pix - 1) + totalSamps)%totalSamps];
            var dif1re = timeShifts[minIndex*totalSamps + (2*pix + totalSamps)%totalSamps];
            var dif2re = timeShifts[minIndex*totalSamps + (2*(pix + 1) + totalSamps)%totalSamps];
            */
            /*
            var dif0im = timeShifts[minIndex*totalSamps + (2*(pix - 1) + totalSamps +1 )%totalSamps];
            var dif1im = timeShifts[minIndex*totalSamps + (2*pix + totalSamps +1)%totalSamps];
            var dif2im = timeShifts[minIndex*totalSamps + (2*(pix + 1) + totalSamps +1)%totalSamps];
            */
            var sh = shifts[minIndex]/(2*nsh);
            var shift_offset = minIndex*totalSamps;
            if(sh>0){
                var dif0re = timeShifts[shift_offset + (2*(pix - 1) + totalSamps)%totalSamps];
                var dif1re = timeShifts[shift_offset + (2*pix + totalSamps)%totalSamps];
                data[line_idx+2*pix] = (dif1re*(1-sh) + dif0re*sh) ;
                //data[line_idx+2*pix +1] = (dif1im*(1-sh) + dif0im*sh) ;
            }else{
                var dif1re = timeShifts[shift_offset + (2*pix + totalSamps)%totalSamps];
                var dif2re = timeShifts[shift_offset + (2*(pix + 1) + totalSamps)%totalSamps];
                //sh= 0-sh;
                data[line_idx+2*pix] = (dif1re*(1+sh) - dif2re*sh) ;
                //data[line_idx+2*pix +1] = (dif1im*(1-sh) + dif2im*sh) ;
            }
            
        }//for pix
        //DEBUG
        //end_time = performance.now();
        //calc_time += end_time-start_time;
    }//for line


    //DEBUG
    //console.log('shift time:',shift_time,'calc time:',calc_time);
}

function exportnii(vec,filename){
    fs.writeFileSync(filename,vec.toString());
}
