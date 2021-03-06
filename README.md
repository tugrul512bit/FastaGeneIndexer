# FastaGeneIndexer
A compressed nucleobase sequence cache for FASTA formatted files, backed by the combined graphics-card-memory of system with 1 microsecond sub-sequence access latency and 300+MB/s throughput using multithreaded queries, on a low-end computer (fx8150 cpu + 3x low-end gpus).

Wiki: https://github.com/tugrul512bit/FastaGeneIndexer/wiki/Sub-Sequence-Query

Dependencies:

https://github.com/tugrul512bit/VirtualMultiArray

C++17

At least one graphics card that supports OpenCL + developer tools (to compile)

At least one graphics card that supports OpenCL + driver (to run)

How it works:

It parses a ".fna" or ".faa" FASTA gene sequence file and indexes all parts into video-memory buffers to get high random-access performance (that scales well with number of graphics cards and cpu threads in use) to all indices. Saving RAM + disk usage lets program stream relevant "result" data to disk without stuttering, unless all video memory is needed elsewhere.

Compression algorithm (not really compression but an encoding that reduces size of text) is Huffman Encoding so it suits well to any type of sequence (and its descriptor). Compression ratio is between 1.3x - 3.7x.

For a system with 3.6GHz FX8150 + 3 low-end graphics cards (2GB memory each), this can store(and access) ~20GB base pair data without using more than 2 GB RAM, at 250-300 MB/s throughput (multithreaded&cached readings) and 4-17microseconds latency (multi-threaded - single threaded), without touching disk drives (but after a "slower" index initialization is complete).

Compiler flags for example main.cpp:

```
g++-10 -std=c++17 -O3 -g3 -Wall -c -fmessage-length=0 -mavx -march=native 
-funroll-loops -mtune=native -funsafe-math-optimizations -ftree-vectorize 
-ffast-math  -fopenmp -pthread -fPIC -MMD -MP -MF"src/main.d" -MT"src/main.d" 
-o "src/main.o" "../src/main.cpp"
```

Linker flags:
```
g++-10 -L/usr/local/cuda/lib64 -o "Virtualizer2"  ./src/main.o   -lOpenCL -lgomp -lpthread
```

How it is used:

```cpp
try
{
	bool debug=true; // to see each pass of encoding, Huffman tree leaf values, compression amount
	FastaGeneIndexer fasta("influenza.fna",debug); // a 1.4GB file that goes down to 428MB and shared between 3 gpus, ~142MB each, (RAM: 250MB used)
	std::string sequence = fasta.getSequence(817586);// ACTG...
	std::cout<<sequence<<std::endl;
}
catch (std::exception& ex)
{
	std::cout << ex.what();
}
```

```
gi|1905469918|gb|MH981935|Influenza A virus (A/mallard/Alberta/273/2017(H4N6)) segment 1 polymerase PB2 (PB2) gene, complete cd6
ATGGAGAGAATAAAAGAACTAAGAGATCTAATGTCACAGTCTCGCACTCGCGAGATACTAACCAAAACCACTGTTGACCACATGGCCATAATCAAAAAGTACACGTCGGGAAGACAAGAGAAGAACCCCGCACTCAGGATGAAATGGATGATGGCAATGAAATATCCAATCACAGCGGATAAGCGAATAATGGAAATGATCCCTGAAAGGAATGAACAAGGGCAAACCCTCTGGAGCAAAACAAACGATGCCGGATCAGACCGAGTGATGGTATCACCTCTGGCTGTGACATGGTGGAATAGGAATGGACCAACAACGAGTACAGTTCATTATCCAAAGGTATATAAAACTTATTTCGAAAAAGTCGAAAGGTTGAAACATGGGACCTTTGGCCCTGTACACTTCAGGAACCAAGTTAAGATAAGACGGAGGGTCGACATAAACCCGGGCCATGCTGACCTCAGTGCCAAAGAGGCACAGGATGTAATCATGGAAGTTGTCTTCCCAAATGAAGTGGGAGCGAGAATACTGACGTCGGAGTCACAACTGACGATAACAAAGGAGAAGAAGGAAGAACTCCAGGACTGCAAAATCGCCCCTCTGATGGTTGCATACATGCTAGAAAGAGAGTTGGTCCGCAAGACGAGGTTTCTCCCAGTGGCTGGTGGAACAAGCAGTGTCTACATTGAAGTGCTGCATTTGACCCAGGGGACATGCTGGGAGCAGATGTACACTCCAGGAGGAGAAGTGAGAAACGATGATGTAGACCAGAGCTTGATCATTGCTGCCAGGAACATAGTAAGAAGAGCAACGGTATCAGCAGACCCACTAGCATCTCTATTGGAGATGTGCCACAGCACACAAATTGGGGGAATAAGAATGGTAGACATTCTTCGGCAAAATCCCACAGAGGAACAAGCCGTGGACATATGCAAGGCAGCAATGGGCTTGAGAATTAGCTCATCTTTCAGTTTCGGTGGATTCACTTTTAAAAGAACAAGTGGATCGTCAGTCAAAAGGGAAGAAGAAGTGCTTACGGGCAACCTTCAAACATTGAAAATAAGAGTACATGAGGGATATGAAGAGTTTACAATGGTCGGAAGAAGAGCAACGGCCATTCTCAGGAAGGCAACCAGAAGACTGATCCAGCTAATAGTAAGTGGAAGAGATGAACAGTCAATTGCTGAAGCAATAATTGTGGCCATGGTATTCTCACAAGAGGACTGCATGATTAAGGCAGTTCGAGGTGATCTGAATTTTGTCAATAGAGCGAACCAGCGGTTGAACCCAATGCATCAGCTCTTGAGACACTTTCAAAAGGATGCAAAAGTGCTTTTCCAAAATTGGGGAATCGAACCCATTGACAATGTGATGGGAATGATCGGGATATTGCCTGACATGACCCCGAGTACTGAGATGTCGCTGAGGGGAATAAGAGTCAGTAAGATGGGAGTAGATGAATACTCCAGCACAGAACGGGTAGTAGTAAGCATCGACCGATTTTTAAGAGTTCGAGATCAACGGGGGAACGTACTATTGTCACCCGAGGAAGTCAGCGAGACACAAGGAACAGAGAAACTGACAATCACTTATTCGTCATCAATGATGTGGGAGATCAATGGCCCTGAGTCGGTGTTGGTCAATACTTATCAGTGGATAATTAGAAACTGGGAAACTGTAAAAATTCAATGGTCACAGGATCCCACAATGCTGTATAATAAGATGGAATTCGAGCCATTTCAGTCTCTGGTCCCTAAGGCAGCCAGAGGTCAATACAGTGGGTTCGTGAGGACACTATTCCAGCAAATGCGAGATGTGCTTGGAACATTTGACACTGTTCAGATAATAAAACTACTTCCCTTTGCTGCTGCCCCACCGGAACAAAGTAGGATGCAGTTCTCCTCTCTGACTGTGAATGTAAGAGGATCAGGAATGAGAATACTAGTAAGAGGCAATTCCCCAGTGTTCAATTACAACAAGGCCACCAAGAGGCTCACAGTTCTCGGGAAAGATGCAGGTGCATTGACAGAAGATCCAGATGAAGGCACAGCTGGGGTGGAGTCCGCTGTTTTAAGAGGATTCCTCATTTTGGGCAAAGAAGACAAGAGATATGGCCCAGCATTGAGCATCAATGAGCTGAGCAATCTTGCAAAGGGAGAGAAGGCTAATGTGCTAATTGGGCAAGGAGACGTGGTGTTGGTGATGAAACGGAAACGGGACTCTAGCATACTTACTGACAGCCAGACAGCGACCAAAAGAATTCGGATGGCCATCAATTAG

```

Benchmark:

```cpp
#include "FastaGeneIndexer.h"
#include "lib/CpuBenchmarker.h"
int main(int argC, char** argV)
{
    try
    {
        bool debug = true;
        FastaGeneIndexer cache("./data/influenza.fna", debug);

        {
            CpuBenchmarker bench(0,"read 1");
            std::cout<<cache.getSequence(0)<<std::endl;
        }

        {
            CpuBenchmarker bench(0,"read 2 same location");
            cache.getSequence(0);
        }

        // heating cpu
         for(int i=0;i<10000;i++)
         {
             cache.getSequence(i);
         }

         size_t count=0;
         {
             CpuBenchmarker bench(0,"read 10k sequences, single thread",10000);
             for(int i=0;i<10000;i++)
             {
                 count += cache.getSequence(i).size();
             }
         }

         {
             CpuBenchmarker bench(count,"read 10k sequences, multi-thread",10000);

             #pragma omp parallel for
             for(int i=0;i<10000;i++)
             {
                 cache.getSequence(i);
             }
         }
	 
	 {
             CpuBenchmarker bench(50000,"read 10k sub-sequences(5 symbols each), single-thread",10000);

             for(int i=0;i<10000;i++)
             {
                 cache.getSequence(i,5,5);
             }
         }
    }
    catch(std::exception & e)
    {
        std::cout<< e.what() <<std::endl;
    }
    return 0;
}

```

Output for FX8150 at 3.6GHz(no turbo) + 1 channel 1333MHz DDR3 RAM (4GB) + 3 low-end 2GB graphics cards (+ test data of https://ftp.ncbi.nih.gov/genomes/INFLUENZA/influenza.fna):

```
AGCAAAAGCAGGAGTTCAAAATGAATCCAAATCAGAAAATAATAACCATTGGGTCAATCTGTATGGGGATCGGAATAATCAGCCTAATATTACAAATTGGAAACATAATCTCAATGTGGGTTAGTCATTCAATTCAGACTGAAAATCAAAATCACCATGAAGCATGCAACCCAAGCATTGCTGGACAGGATGCAGCTTCAGTGGCACTAGCAGGCAATTCCTCTCTTTGTCCCATTAGTGGGTGGGCTATATACAGTAAAGACAATGGTATAAGAATTGGCTCCAAAGGAGACGTATTTGTCATAAGAGAGCCATTTATTTCATGCTCTCACTTGGAATGCAGGACCTTTTTTCTGACTCAGGGCGCCTTGTTGAATGACAAACACTCAAATGGAACCGTTAAAGACAGAAGCCCTTATAGAACCTTGATGAGCTGTCCTGTTGGTGAAGCTCCTTCTCCCTACAATTCAAGGTTCGTGTCGGTTGCATGGTCAGCAAGTGCCTGCCATGATGGCATGGGTTGGCTAACAATCGGAATTTCTGGTCCAGATAATGGAGCAGTGGCTGTATTAAAATACAATGGTATAATAACAGACACCATCAAAAGTTGGAAAAATAACATATTGAGAACGCAAGAGTCTGAATGTGCCTGTATAAATGGTTCATGTTTCACTATAATGACCGATGGCCCAAGTAATGGGCAGGCCTCGTACAAAATTTTCAAGATAGAAAAGGGGAAGGTAGTCAAATCAAGTGAATTAAATGCACCTAATTACCACTACGAGGAATGTTCTTGTTATCCTGATGCAGGTGAAGTAATGTGTGTATGCAGGGATAATTGGCATGGTTCGAATCGACCATGGGTGTCTTTCAATAAAAACCTTGATTATCAAATAGGTTACATCTGCAGTGGGGTTTTTGGTGATAATCCACGACCCAATGATGGAACAGGCAGCTGTGGTCCAGTGTCTTCTAATGGAGCATATGGGATAAAGGGGTTCTCATTTAAGTATGGTAATGGTGTTTGGATAGGGAGAACCAAAAGCACTAGTTCCAGAAGCGGGTTTGAGATGATTTGGGATCCCAATGGATGGACAGAAACTGACAGCAGTTTCTCTGTGAAGCAAGATATTGTAGCAATAACTGATTGGTCGGGATATAGCGGGAGTTTCGTCCAACACCCTGAATTAACAGGGTTGGACTGCATGAGGCCTTGTTTCTGGGTTGAACTGATCAGGGGACGGCCCAACCACAATACGATCTGGACTAGTGGGAGCAGCATTTCCTTCTGTGGTGTGAACAGCGATACTGTAGGTTGGTCTTGGCCAGACGGTGCTGAGTTGCCATTCACCATTGACAAGTAGTTTGTTCAAAAAAACTCCTTGTTTCTACT
read 1: 85321 nanoseconds    
read 2 same location: 26316 nanoseconds    
read 10k sequences, single thread: 176683019 nanoseconds     
(throughput = 17668.30 nanoseconds per iteration) 

read 10k sequences, multi-thread: 37886537 nanoseconds     
(bandwidth = 299.84 MB/s)      (throughput = 3788.65 nanoseconds per iteration) 

read 10k sub-sequences(5 symbols each), single-thread: 11770656 nanoseconds     
(bandwidth = 4.25 MB/s)      (throughput = 1177.07 nanoseconds per iteration)

```

Just reading 5 nucleobase symbols is as fast as 1million queries per second per CPU thread, on a slow development machine. If queries are in page cache (in the internal virtual array), then latency drops to ~100 nanoseconds (also multi-threaded bandwidth increases to ~400 MB/s).
