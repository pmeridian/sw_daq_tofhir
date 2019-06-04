#ifndef __PETSYS__EVENT_DECODE_HPP__DEFINED__
#define __PETSYS__EVENT_DECODE_HPP__DEFINED__

#include <stdint.h>
#include <string>
#include <stdio.h>

namespace PETSYS {       

class RawEventWord{

public:
	RawEventWord(unsigned __int128 word) : word(word){};
	~RawEventWord() {};
	
	unsigned short getQFine()       { return (word >> 0) % 1024; };
	unsigned short getQCoarse()     { return (word >> 10) % 1024; };
	unsigned short getT2Fine()      { return (word >> 20) % 1024; };
	unsigned short getT2Coarse()    { return (word >> 30) % 1024; };
	unsigned short getT1Fine()      { return (word >> 40) % 1024; };
	unsigned short getT1Coarse()    { return (word >> 50) % 1024; };
	unsigned short getIdleTime()    { return (word >> 60) % 1024; };
	unsigned short getTriggerBits() { return (word >> 70) % 16; };
	unsigned short getTacID()      { return (word >> 74) % 4; };
	unsigned int getChannelID()   { return word >> 76; };
	
	
// 	unsigned getEFine() {
// 		unsigned v = (word >> 20) % 1024;
// 		return v;
// 	};
// 
// 	unsigned getTFine() {
// 		unsigned v = (word>>40) % 1024;
// 		return v;
// 	};
// 
// 	unsigned getECoarse() {
// 		return (word>>30) % 1024;
// 	};
// 	
// 	unsigned getTCoarse() {
// 		return (word>>50) % 1024;
// 	};
// 	
// 	unsigned getTacID() {
// 	
// 		return (word>>74) % 4;
// 	};
// 	
// 	unsigned getChannelID() {
// 	
// 		return word>>76;
// 	};
	
public:
	unsigned __int128 word;

};

}
#endif
