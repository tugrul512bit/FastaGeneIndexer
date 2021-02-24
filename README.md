# FastaGeneIndexer
A compressed FASTA sequence cache backed by the combined video memory of system.

Dependencies:

https://github.com/tugrul512bit/VirtualMultiArray

C++17

At least one graphics card that supports OpenCL + developer tools (to compile)

At least one graphics card that supports OpenCL + driver (to run)

How it works:

It parses a ".fna" or ".faa" FASTA gene sequence file and indexes all parts into video-memory buffers to get high random-access performance (that scales well with number of graphics cards and cpu threads in use) to all indices. Saving RAM + disk usage lets program stream relevant "result" data to disk without stuttering, unless all video memory is needed elsewhere.

Compression algorithm is Huffman Encoding so it suits well to any type of sequence (and its descriptor).
