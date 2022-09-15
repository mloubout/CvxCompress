# CvxCompres library
This CvxCompress was originally authored by Thor Johnsen.    

# Building
## GCC

It is believed that this will work with a wide variety of gcc
versions, but as I am writing this I have only tested it with gcc 7.30 and gcc 9.3.0. 

## Tests

There are several tests included, some of which depends on binary
cubes of input need to exist at some absolute path which is no longer
present. Thus, here I describe only the ones that require no input file. 

```
Prompt> ( setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${cwd} ; ./CvxCompress_Test )
*
* CvxCompress module tests (ICC1900, AVX).
*

2. Verify correctness of forward wavelet transform...[Passed]
3. Verify correctness of inverse wavelet transform...[Passed]
4. Test throughput of wavelet transform (forward + inverse)...
 ->   8 x   8 x   8 ( L1 ) ::  4.770 secs - 3602 MCells/s - 435 GF/s
 ->  16 x  16 x  16 ( L1 ) ::  4.629 secs - 3711 MCells/s - 480 GF/s
 ->  32 x  32 x  32 ( L2 ) ::  6.371 secs - 2697 MCells/s - 361 GF/s
 ->  64 x  64 x  64 ( L3 ) ::  7.438 secs - 2310 MCells/s - 314 GF/s
 -> 128 x 128 x 128 (DRAM) :: 11.056 secs - 1557 MCells/s - 213 GF/s
 -> 256 x 256 x 256 (DRAM) :: 18.059 secs - 965 MCells/s - 133 GF/s

5. Verify correctness of Copy_To_Block method...[Passed!]
6. Verify correctness of Copy_From_Block method...[Passed!]
7. Test throughput of block copy...
 ->   8 x   8 x   8 [Passed!] ::  0.413 secs - 2603 MCells/s - 31.23 GB/s
 ->  16 x  16 x  16 [Passed!] ::  0.438 secs - 2454 MCells/s - 29.44 GB/s
 ->  32 x  32 x  32 [Passed!] ::  0.605 secs - 1775 MCells/s - 21.30 GB/s
 ->  64 x  64 x  64 [Passed!] ::  0.517 secs - 2077 MCells/s - 24.92 GB/s
 -> 128 x 128 x 128 [Passed!] ::  0.608 secs - 1767 MCells/s - 21.20 GB/s
 -> 256 x 256 x 256 [Passed!] ::  0.694 secs - 1547 MCells/s - 18.56 GB/s

8. Verify correctness of Global_RMS method...[Passed!]
9. Test throughput of Compress() method...
nx=1604, ny=1606, nz=968
 ->   8 x   8 x   8 ( L1 )  9 iterations - 11.024 secs - 2036 MCells/s - ratio 29.03:1
 ->  16 x  16 x  16 ( L1 ) 12 iterations - 10.419 secs - 2872 MCells/s - ratio 51.37:1
 ->  32 x  32 x  32 ( L2 ) 10 iterations - 10.356 secs - 2408 MCells/s - ratio 73.40:1
 ->  64 x  64 x  64 ( L3 ) 10 iterations - 10.654 secs - 2341 MCells/s - ratio 93.21:1
 -> 128 x 128 x 128 (DRAM)  7 iterations - 10.751 secs - 1624 MCells/s - ratio 105.98:1
10. Test throughput of Decompress() method...
 ->   8 x   8 x   8 ( L1 )  7 iterations - 10.321 secs - 1691 MCells/s
 ->  16 x  16 x  16 ( L1 ) 10 iterations - 10.283 secs - 2425 MCells/s
 ->  32 x  32 x  32 ( L2 )  8 iterations - 11.327 secs - 1761 MCells/s
 ->  64 x  64 x  64 ( L3 )  8 iterations - 10.695 secs - 1865 MCells/s

Prompt> ( setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:$cwd ; ./Test_With_Generated_Input )
Using global RMS.
Using size of input (nslow, middle, fast) = ( 320, 416, 352 ) 
Starting compression test
RMS:
 input :     0.707107
 output:     0.707080
 Difference: 0.000118
compression ratio (return value) = 1048.35:1, compression throughput = 1282 MC/s, decompression throughput = 1207 MC/s, error = 1.663531e-04, SNR = 75.6 dB
Total compression and decompression times were 0.08 seconds
Total compression ratio (based on compress_length) was 1048.35:1, compressed length in bytes = 178789 
Using size of input (nslow, middle, fast) = ( 640, 832, 704 ) 
Starting compression test
RMS:
 input :     0.707107
 output:     0.707080
 Difference: 0.000118
compression ratio (return value) = 1048.55:1, compression throughput = 1838 MC/s, decompression throughput = 1619 MC/s, error = 1.663531e-04, SNR = 75.6 dB
Total compression and decompression times were 0.44 seconds
Total compression ratio (based on compress_length) was 1048.55:1, compressed length in bytes = 1430039 
Using size of input (nslow, middle, fast) = ( 960, 1248, 1056 ) 
Starting compression test
RMS:
 input :     0.707107
 output:     0.707080
 Difference: 0.000118
compression ratio (return value) = 1048.57:1, compression throughput = 2217 MC/s, decompression throughput = 1708 MC/s, error = 1.663531e-04, SNR = 75.6 dB
```

