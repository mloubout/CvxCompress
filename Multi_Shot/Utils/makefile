all: MergeShotAndReceiverASCIIFiles DumpArrND ExtractReceiver DumpTraceHeaders_Custom1 DayOfYearToDate Interactive_Coordinate_Converter_rev1
	cp ./MergeShotAndReceiverASCIIFiles ../../bin
	cp ./DumpArrND ../../bin
	cp ./ExtractReceiver ../../bin
	cp ./DumpTraceHeaders_Custom1 ../../bin
	cp ./DayOfYearToDate ../../bin
	cp ./Interactive_Coordinate_Converter_rev1 ../../bin

MergeShotAndReceiverASCIIFiles: MergeShotAndReceiverASCIIFiles.cpp
	icc -std=c++11 -I../Geom -I../ArrND -I../../ $^ -L../../lib -lvoxet -o $@

ExtractReceiver: ExtractReceiver.cpp
	icc -std=c++11 -I../Geom -I../ArrND -I../../Common $^ -L../../lib -lvoxet -o $@

DumpArrND: DumpArrND.cpp
	icc -std=c++11 -I../Geom -I../ArrND -I../../Common $^ -L../../lib -lvoxet -o $@

DumpTraceHeaders_Custom1: DumpTraceHeaders_Custom1.cpp
	icc -std=c++11 -I../../Common $^ -L../../lib -lsegy -lcrypt -o $@

DayOfYearToDate: DayOfYearToDate.cpp
	icc -std=c++11 $^ -o $@

Interactive_Coordinate_Converter_rev1: Interactive_Coordinate_Converter_rev1.cpp
	icc -std=c++11 -I../Geom -I../ArrND -I../../Common $^ -L../../lib -lvoxet -o $@
