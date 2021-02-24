# FastaGeneIndexer
A compressed FASTA sequence cache backed by the combined video memory of system.

Dependencies:

https://github.com/tugrul512bit/VirtualMultiArray

C++17

At least one graphics card that supports OpenCL + developer tools (to compile)

At least one graphics card that supports OpenCL + driver (to run)

How it works:

It parses a ".fna" or ".faa" FASTA gene sequence file and indexes all parts into video-memory buffers to get high random-access performance (that scales well with number of graphics cards and cpu threads in use) to all indices. Saving RAM + disk usage lets program stream relevant "result" data to disk without stuttering, unless all video memory is needed elsewhere.

Compression algorithm is Huffman Encoding so it suits well to any type of sequence (and its descriptor). Compression ratio is between 1.3x - 3.7x.

For a system with 2.1GHz FX8150 + 3 low-end graphics cards (2GB memory each), this can store(and access) ~20GB base pair data without using more than 1 GB RAM, at 250-300 MB/s throughput (multithreaded&cached readings) and 15microseconds latency (single threaded un-cached access), without touching disk drives (after indexing is complete).
