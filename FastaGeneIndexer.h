/*
 * FastaGeneIndexer.h
 *
 *  Created on: Feb 21, 2021
 *      Author: tugrul
 */

#ifndef FASTAGENEINDEXER_H_
#define FASTAGENEINDEXER_H_

#include "../GraphicsCardSupplyDepot.h"
#include "../VirtualMultiArray.h"
#include <fstream>
#include <algorithm>
#include <omp.h>
#include <fstream>
#include <string>
#include <filesystem>
#include <mutex>
#include <vector>
#include <map>


// C++17


inline void storeBit(unsigned char & data, const bool value, const int pos) noexcept
{
	data = (value << pos) | (data & ~(1 << pos));
}

inline bool loadBit(const unsigned char & data, const int pos) noexcept
{
	return (data>>pos)&1;
}


struct Node
{
	size_t count;
	int self;
	int leaf1;
	int leaf2;
	unsigned char data;
	bool isLeaf;
};


class HuffmanTree
{


public:
	HuffmanTree(){
		for(int i=0;i<256;i++)
		{
			referenceMapDirect[i]=0;
			encodeDirectSize[i]=0;
		}

	}

	void add(unsigned char data)
	{
		referenceMapDirect[data]++;

	}

	void generateTree(const bool debug=false)
	{
		std::vector<Node> sortedNodes;


		int ctr=0;

		for(int i=0;i<256;i++)
		{
			size_t ct = referenceMapDirect[i];
			if(ct>0)
			{
				Node node;

				node.data=i;
				node.count=ct;

				node.self=ctr;
				node.leaf1=-1;
				node.leaf2=-1;
				node.isLeaf=true;
				referenceVec.push_back(node);
				sortedNodes.push_back(node);
				ctr++;
			}
		}

		std::sort(sortedNodes.begin(), sortedNodes.end(),[](const Node & n1, const Node & n2){ return n1.count<n2.count;});

		while(sortedNodes.size()>1)
		{
			Node node1 = sortedNodes[0];
			Node node2 = sortedNodes[1];
			Node newNode;
			newNode.count = node1.count + node2.count;
			newNode.data=0;
			newNode.leaf1 = node1.self;
			newNode.leaf2 = node2.self;
			newNode.self = ctr;
			newNode.isLeaf=false;
			sortedNodes.erase(sortedNodes.begin());
			sortedNodes.erase(sortedNodes.begin());
			sortedNodes.push_back(newNode);

			referenceVec.push_back(newNode);
			std::sort(sortedNodes.begin(), sortedNodes.end(),[](const Node & n1, const Node & n2){ return n1.count<n2.count;});
			ctr++;
		}

		root = sortedNodes[0];


		std::function<void(Node,std::vector<bool>)> g = [&](Node node, std::vector<bool> path){
			if(node.leaf1!=-1)
			{
				std::vector<bool> path1 = path;
				path1.push_back(false);
				g(referenceVec[node.leaf1],path1);
			}

			if(node.leaf2!=-1)
			{
				std::vector<bool> path2 = path;
				path2.push_back(true);
				g(referenceVec[node.leaf2],path2);
			}

			if((node.leaf1 == -1) && (node.leaf2 == -1))
			{
				encodeDirect[node.data]=path;
				encodeDirectSize[node.data]=path.size();
			}

		};

		std::vector<bool> rootPath;
		g(root,rootPath);

		if(debug)
		{
			std::cout<<"-------------------------------------"<<std::endl;

			for(int i=0;i<256;i++)
			{
				if(encodeDirect[i].size()>0)
				{
					std::cout<<(unsigned char)i<<": ";
					for(const auto & f: encodeDirect[i])
					{
						std::cout<<f<<" ";
					}
					std::cout<<std::endl;
				}
			}
			std::cout<<"-------------------------------------"<<std::endl;
		}
	}

	size_t getCount(unsigned char data)
	{
		return referenceMapDirect[data];
	}

	inline
	const std::vector<bool> & generateBits(unsigned char data) const noexcept
	{
		return encodeDirect[data];
	}

	inline
	const int & generateBitsSize(unsigned char data) const noexcept
	{
		return encodeDirectSize[data];
	}


	inline
	const unsigned char followBitsDirect(const Node * __restrict__ const refVec, const unsigned char * __restrict__ const path, size_t & idx, const size_t & ofs) const
	{
		unsigned char result;
		const Node * curNode=&root;
		bool work=true;

		while(work)
		{
			int p = loadBit(path[(idx>>3)-ofs],idx&7);
			if(curNode->isLeaf)
			{
				result=curNode->data;
				work=false;
			}
			else
			{
				curNode = refVec+(curNode->leaf2*p + curNode->leaf1*(1-p));
				idx++;
			}
		}
		return result;
	}

	const Node * getRefData() const { return referenceVec.data(); }

	~HuffmanTree(){}
private:
	Node root;
	std::vector<Node> referenceVec;
	std::map<unsigned char,size_t> referenceMap;
	size_t referenceMapDirect[256];
	std::vector<bool> encodeDirect[256];
	int encodeDirectSize[256];
};

// FASTA file indexer that caches bits of data in video-memory
// supports maximum 12 physical graphics cards (with combined video memory size of them)
// supports 10 million indices
// uses RAM temporarily = (total video memory size in GB) * 80MB = 100 GB total VRAM means 8GB temporary RAM usage
class FastaGeneIndexer
{
public:
	FastaGeneIndexer(){}
	FastaGeneIndexer(std::string fileName, bool debug=false){
		fileDescriptorN=0;
		sizeIO = 1024*1024*16;
		// get file size
		size_t bytes = countFileBytes(fileName);

		// first-pass for Huffman encoding: generate Huffman tree
		// second-pass for Huffman encoding: count total bits for encoded file
		// third-pass: fill virtual array with encoded bits


		if(debug)
		{
			std::cout<<"total bytes (padded to 16-MB multiple): "<<bytes<<std::endl;
			std::cout<<"Generating Huffman Tree for descriptors and sequences."<<std::endl;
		}

		// Huffman Tree generation
		generateHuffmanTree(fileName,bytes,debug);

		if(debug)
			std::cout<<"Counting number of bits to produce."<<std::endl;

		// bit counting for processed file
		size_t totalBits = countBits(fileName,bytes,debug);

		if(debug)
		{
			std::cout<<"descriptor: "<<descriptorBits<<" bits, sequence: "<<sequenceBits<<" bits"<<std::endl;
			std::cout<<"Compression: from "<<(bytes/1000000.0)<<" MB (padded to 16MB-multiple) to "<<(totalBits/8 + 1)/1000000.0<<" MB"<<std::endl;
			std::cout<<"Allocating virtual array."<<std::endl;
		}



		// alllocate
		size_t totalBytes = (totalBits/8 + 1);


		{
			size_t bytesPerSequence = (totalBytes/fileDescriptorN);
			size_t bytesNeededPerPage = 10*bytesPerSequence;
			pageSize=128;
			while(pageSize<bytesNeededPerPage)
			{
				pageSize *= 2;
			}

			while(pageSize > sizeIO)
			{
				pageSize /= 2;
			}

			while(pageSize > (totalBytes / 12))
			{
				pageSize /= 2;
			}

			int numCachePerGpu = sizeIO / pageSize;
			if(numCachePerGpu>15)
				numCachePerGpu=15;


			size_t n = totalBytes + pageSize - (totalBytes%pageSize);

			if(debug)
			{
				std::cout<<"page size for virtual array = "<<pageSize<<" bytes"<<std::endl;
				std::cout<<"virtual array size = "<<n<<" bytes"<<std::endl;
				std::cout<<"file i/o size = "<<sizeIO<<" bytes"<<std::endl;
			}
			GraphicsCardSupplyDepot gpu;
			data = VirtualMultiArray<unsigned char>(n,gpu.requestGpus(),pageSize,numCachePerGpu,{1,1,1,1,1,1,1,1,1,1,1,1},VirtualMultiArray<unsigned char>::MemMult::UseVramRatios);
		}

		if(debug)
			std::cout<<"Filling virtual array with bits of file."<<std::endl;

		{
			descriptorBeginBit.resize(fileDescriptorN);
			descriptorBitLength.resize(fileDescriptorN);
			sequenceBeginBit.resize(fileDescriptorN);
			sequenceBitLength.resize(fileDescriptorN);

			size_t writtenBytes = readFileIntoArray(fileName, bytes, debug);

			if(debug)
			{
				std::cout<<"Written "<<writtenBytes<<" bytes ("<<writtenBytes/1000000.0<<" MB) to virtual array"<<std::endl;
				std::cout<<"Descriptor: "<<descriptorBeginBit.size()<<" "<<descriptorBitLength.size()<<std::endl;
				std::cout<<"Sequence: "<<sequenceBeginBit.size()<<" "<<sequenceBitLength.size()<<std::endl;
			}
		}
	}

	// get a gene descriptor at index=id without line-feed nor '>' characters
	// thread-safe
	std::string getDescriptor(size_t id)
	{
		std::string result;
		size_t i0 = descriptorBeginBit[id];
		unsigned int r0 = descriptorBitLength[id];
		const size_t i03 = (i0>>3);
		const int r1 = (r0>>3) + 1;
		std::vector<unsigned char> tmpData = data.readOnlyGetN(i03,r1);
		size_t pos = i0;
		const unsigned int pL = r0+pos;
		const unsigned char * dt = tmpData.data();
		const Node * nodePtr = descriptorCompression.getRefData();
		while(pos<pL)
		{
			result += descriptorCompression.followBitsDirect(nodePtr,dt,pos,i03);
		}
		return result;
	}

	// get a gene sequence at index=id without line-feed characters
	// thread-safe
	std::string getSequence(size_t id)
	{
		std::string result;
		size_t i0 = sequenceBeginBit[id];
		unsigned int r0 = sequenceBitLength[id];
		const size_t i03 = (i0>>3);
		const int r1 = (r0>>3) + 1;
		std::vector<unsigned char> tmpData = data.readOnlyGetN(i03,r1);
		size_t pos = i0;
		const unsigned int pL = r0+pos;
		const unsigned char * dt = tmpData.data();
		const Node * nodePtr = sequenceCompression.getRefData();
		while(pos<pL)
		{
			result += sequenceCompression.followBitsDirect(nodePtr,dt,pos,i03);
		}
		return result;
	}

	// this is required for getSequenceByDescriptor()
	// maps descriptor strings to index values
	void initDescriptorIndexMapping()
	{
		if(descriptorToIndex.size()==0)
		{
			const size_t nSeq = n();

			for(size_t i=0;i<nSeq;i++)
			{
				descriptorToIndex[getDescriptor(i)]=i;
			}
		}
	}

	// get a gene sequence by its descriptor string
	// thread-safe
	std::string getSequenceByDescriptor(std::string name)
	{
		return getSequence(descriptorToIndex[name]);
	}

	// number of sequences
	size_t n()
	{
		return descriptorBeginBit.size();
	}

	// returns total number of occurences of a symbol 'A', 'C', 'G', ...
	size_t getSymbolCount(unsigned char symbol)
	{
		return sequenceCompression.getCount(symbol);
	}
private:

	// byte data in video memories
	VirtualMultiArray<unsigned char> data;
	HuffmanTree descriptorCompression;
	HuffmanTree sequenceCompression;
	size_t descriptorBits;
	size_t sequenceBits;
	std::vector<size_t> descriptorBeginBit;
	std::vector<size_t> sequenceBeginBit;
	std::vector<int> descriptorBitLength;
	std::vector<int> sequenceBitLength;
	size_t fileDescriptorN;
	size_t pageSize;
	size_t sizeIO;
	std::map<std::string,size_t> descriptorToIndex;

	// returns file size in resolution of 1024*1024*16 bytes (for the paging performance of virtual array)
	// will require to set '\0' for excessive bytes of last block
	size_t countFileBytes(std::string inFile)
	{
		size_t size = std::filesystem::file_size(inFile);
		return size + sizeIO - ( size%(sizeIO) );
	}


	void generateHuffmanTree(std::string inFile, size_t bytes, const bool debug=false)
	{
		std::ifstream bigFile(inFile);
		const size_t bufferSize = sizeIO;
		std::vector<unsigned char> buf;
		buf.resize(bufferSize);
		size_t ctr=0;
		const int div = bytes/bufferSize;
		int ctrDebug = 0;
		bool encodingDescriptor = false;
		bool encodingSequence = false;
		int encodedLength = 0;
		while(bigFile)
		{
			if(debug)
			{

				if(ctrDebug==1)
				{
					ctrDebug=0;
					std::cout<<"Progress: "<<(ctr/bufferSize)<<"/"<<div<<std::endl;
				}
				ctrDebug++;
			}

			// nullify overflowed bytes
			if(ctr>bytes - bufferSize*2)
			{
				for(int i=0;i<bufferSize;i++)
					buf[i]='\0';
			}

			bigFile.read((char *)buf.data(), bufferSize);



			for(int i=0;i<bufferSize;i++)
			{

				const unsigned char elm = buf[i];

				// if descriptor sign found
				if((elm=='>')/* && (!encodingDescriptor)*/)
				{
					fileDescriptorN++;

					// enable descriptor tree building
					encodingDescriptor = true;
					encodingSequence = false;
					encodedLength = 0;
					continue;
				}

				if((elm=='\n' || elm=='\r' || elm=='\0') && encodingDescriptor && encodedLength>0)
				{
					encodingDescriptor = false;
					encodingSequence = true;
					encodedLength = 0;
					continue;
				}

				if((elm=='\0') && encodingSequence && encodedLength>0)
				{
					encodingDescriptor = false;
					encodingSequence = false;
					encodedLength = 0;
					continue;
				}

				if(elm=='\n' || elm=='\r')
				{
					continue;
				}

				if(encodingDescriptor)
				{
					descriptorCompression.add(elm);
				}

				if(encodingSequence)
				{
					sequenceCompression.add(elm);
				}

				encodedLength++;
			}

			ctr+=bufferSize;
		}

		descriptorCompression.generateTree(debug);
		sequenceCompression.generateTree(debug);
	}

	size_t countBits(std::string inFile, size_t bytes, const bool debug=false)
	{
		sequenceBits=0;
		descriptorBits=0;
		std::ifstream bigFile(inFile);
		const size_t bufferSize = sizeIO;
		std::vector<unsigned char> buf;
		buf.resize(bufferSize);
		size_t ctr=0;

		const int div = bytes/bufferSize;
		int ctrDebug = 0;
		bool encodingDescriptor = false;
		bool encodingSequence = false;
		int encodedLength = 0;
		while(bigFile)
		{
			if(debug)
			{

				if(ctrDebug==1)
				{
					ctrDebug=0;
					std::cout<<"Progress: "<<(ctr/bufferSize)<<"/"<<div<<std::endl;
				}
				ctrDebug++;
			}

			// nullify overflowed bytes
			if(ctr>bytes - bufferSize*2)
			{
				for(int i=0;i<bufferSize;i++)
					buf[i]='\0';
			}

			bigFile.read((char *)buf.data(), bufferSize);


			for(int i=0;i<bufferSize;i++)
			{
				const unsigned char elm = buf[i];

				// if descriptor sign found
				if((elm=='>')/* && (!encodingDescriptor)*/)
				{
					// enable descriptor tree building
					encodingDescriptor = true;
					encodingSequence = false;
					encodedLength = 0;
					continue;
				}

				if((elm=='\n' || elm=='\r' || elm=='\0') && encodingDescriptor && encodedLength>0)
				{
					encodingDescriptor = false;
					encodingSequence = true;
					encodedLength = 0;
					continue;
				}

				if((elm=='\0') && encodingSequence && encodedLength>0)
				{
					encodingDescriptor = false;
					encodingSequence = false;
					encodedLength = 0;
					continue;
				}

				if(elm=='\n' || elm=='\r')
				{
					continue;
				}

				if(encodingDescriptor)
				{
					descriptorBits += descriptorCompression.generateBitsSize(elm);

				}

				if(encodingSequence)
				{
					sequenceBits += sequenceCompression.generateBitsSize(elm);

				}

				encodedLength++;


			}

			ctr+=bufferSize;
		}
		return descriptorBits + sequenceBits;
	}


	size_t readFileIntoArray(std::string inFile, size_t bytes, const bool debug=false)
	{
		size_t currentBit = 0;
		size_t writtenBytes = 0;
		std::ifstream bigFile(inFile);
		const size_t bufferSize = sizeIO;
		std::vector<unsigned char> buf;
		buf.resize(bufferSize);
		size_t ctr=0;


		const int div = bytes/bufferSize;
		int ctrDebug = 0;
		std::vector<unsigned char> encoded;
		encoded.resize(bufferSize);
		int encodedCtr=0;
		unsigned char byteValue=0;
		int byteBitIndex=0;
		bool needsWrite=false;

		bool encodingDescriptor = false;
		bool encodingSequence = false;
		int encodedLength = 0;
		int encodedBitLength = 0;
		size_t descriptorBeginCtr = 0;
		size_t descriptorLengthCtr = 0;
		size_t sequenceBeginCtr = 0;
		size_t sequenceLengthCtr = 0;
		while(bigFile)
		{
			if(debug)
			{

				if(ctrDebug==1)
				{
					ctrDebug=0;
					std::cout<<"Progress: "<<(ctr/bufferSize)<<"/"<<div<<std::endl;
				}
				ctrDebug++;
			}

			// nullify overflowed bytes
			if(ctr>bytes - bufferSize*2)
			{
				for(int i=0;i<bufferSize;i++)
					buf[i]='\0';
			}

			bigFile.read((char *)buf.data(), bufferSize);


			for(int i=0;i<bufferSize;i++)
			{
				const unsigned char elm = buf[i];
				// if descriptor sign found
				if((elm=='>')/* && (!encodingDescriptor)*/)
				{

					// enable descriptor tree building
					if(encodingSequence)
					{
						sequenceBitLength[sequenceLengthCtr++]=encodedBitLength;
					}
					encodedBitLength=0;
					encodingDescriptor = true;
					encodingSequence = false;
					encodedLength = 0;

					descriptorBeginBit[descriptorBeginCtr++]=currentBit;
					continue;
				}

				if(( (elm=='\n') || (elm=='\r') || (elm=='\0')) && encodingDescriptor && (encodedLength>0))
				{
					descriptorBitLength[descriptorLengthCtr++]=encodedBitLength;
					encodedBitLength=0;
					encodingDescriptor = false;
					encodingSequence = true;
					encodedLength = 0;
					sequenceBeginBit[sequenceBeginCtr++]=currentBit;
					continue;
				}

				if((elm=='\0') && encodingSequence && encodedLength>0)
				{
					sequenceBitLength[sequenceLengthCtr++]=encodedBitLength;
					encodedBitLength=0;
					encodingDescriptor = false;
					encodingSequence = false;
					encodedLength = 0;
					continue;
				}

				if(elm=='\n' || elm=='\r')
				{
					continue;
				}

				if(encodingDescriptor)
				{

					const std::vector<bool> & path =descriptorCompression.generateBits(elm);
					const int pL = descriptorCompression.generateBitsSize(elm);
					for(int j=0;j<pL;j++)
					{
						storeBit(byteValue,path[j],byteBitIndex);
						byteBitIndex++;
						needsWrite=true;
						if(byteBitIndex==8)
						{
							// write to byte buffer
							encoded[encodedCtr++]=byteValue;
							byteValue=0;
							byteBitIndex=0;
							needsWrite=false;
						}
					}

					encodedBitLength+=pL;
					currentBit+=pL;
				}

				if(encodingSequence)
				{

					const std::vector<bool> & path = sequenceCompression.generateBits(elm);
					const int pL = sequenceCompression.generateBitsSize(elm);

					for(int j=0;j<pL;j++)
					{
						storeBit(byteValue,path[j],byteBitIndex);
						byteBitIndex++;
						needsWrite=true;
						if(byteBitIndex==8)
						{
							// write to byte buffer
							encoded[encodedCtr++]=byteValue;
							byteValue=0;
							byteBitIndex=0;
							needsWrite=false;
						}
					}

					encodedBitLength+=pL;
					currentBit+=pL;
				}

				encodedLength++;


			}

			size_t eSize = encodedCtr;
			const int nChunks = encodedCtr/pageSize;
			#pragma omp parallel for
			for(int i=0;i<nChunks;i++)
			{
				data.writeOnlySetN(writtenBytes+(i*pageSize),encoded,i*pageSize,pageSize);
			}
			//data.writeOnlySetN(writtenBytes,encoded,0,encodedCtr);
			data.writeOnlySetN(writtenBytes+nChunks*pageSize,encoded,nChunks*pageSize,encodedCtr-nChunks*pageSize);
			writtenBytes += eSize;
			ctr+=bufferSize;
			//encoded.clear();
			encodedCtr=0;
		}

		// don't forget last byte
		if(needsWrite)
		{
			data.set(writtenBytes,byteValue);
			writtenBytes++;
		}

		if(debug)
		{
			std::cout<<descriptorLengthCtr<<" "<<descriptorBeginCtr<<" "<<sequenceLengthCtr<<" "<<sequenceBeginCtr<<std::endl;
		}

		return writtenBytes;
	}
};



#endif /* FASTAGENEINDEXER_H_ */
