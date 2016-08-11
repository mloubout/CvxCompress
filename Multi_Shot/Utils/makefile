all: MergeShotAndReceiverASCIIFiles DumpArrND ExtractReceiver DumpTraceHeaders_Custom1
	cp ./MergeShotAndReceiverASCIIFiles ../../bin
	cp ./DumpArrND ../../bin
	cp ./ExtractReceiver ../../bin
	cp ./DumpTraceHeaders_Custom1 ../../bin
	cp ./DayOfYearToDate ../../bin

MergeShotAndReceiverASCIIFiles: MergeShotAndReceiverASCIIFiles.cpp
	icc -std=c++11 -I../Geom -I../ArrND -I../../ $^ -L../../lib -lvoxet -o $@

ExtractReceiver: ExtractReceiver.cpp
	icc -std=c++11 -I../Geom -I../ArrND -I../../Common $^ -L../../lib -lvoxet -o $@

DumpArrND: DumpArrND.cpp
	icc -std=c++11 -I../Geom -I../ArrND -I../../Common $^ -L../../lib -lvoxet -o $@

DumpTraceHeaders_Custom1: DumpTraceHeaders_Custom1.cpp
	icc -std=c++11 -I../../Common $^ -L../../lib -lsegy -o $@

DayOfYearToDate: DayOfYearToDate.cpp
	icc -std=c++11 $^ -o $@
