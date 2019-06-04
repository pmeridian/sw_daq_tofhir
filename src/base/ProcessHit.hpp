#ifndef __PETSYS__PROCESS_HIT_HPP__DEFINED__
#define __PETSYS__PROCESS_HIT_HPP__DEFINED__

#include <OverlappedEventHandler.hpp>
#include <Event.hpp>
#include <SystemConfig.hpp>
#include <Instrumentation.hpp>


namespace PETSYS {
	
class ProcessHit : public OverlappedEventHandler<RawHit, Hit> {
public:
	static unsigned short remapADC(unsigned short v) { 
		// TOFHiR ADC is bipolar
		// The raw data needs to be remaped to become monotonic
		v = (v < 512) ? (511 - v) : (v);
		
		// TOFHiR ADC has reverse range of TOFPET 2 ADC
		// Artifitially invert it as to re-use calibration code
		v = 1023 - v;
		return v;
	};
	
private:
	SystemConfig *systemConfig;
	bool qdcMode;
	
	uint32_t nReceived;
	uint32_t nReceivedInvalid;
	uint32_t nTDCCalibrationMissing;
	uint32_t nQDCCalibrationMissing;
	uint32_t nXYZMissing;
	uint32_t nSent;
public:
	ProcessHit(SystemConfig *systemConfig, bool qdcMode, EventSink<Hit> *sink);
	virtual void report();
	
protected:
	virtual EventBuffer<Hit> * handleEvents (EventBuffer<RawHit> *inBuffer);
};
	
}

#endif // __PETSYS__PROCESS_HIT_HPP__DEFINED__