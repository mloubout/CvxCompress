all: MergeShotAndReceiverASCIIFiles DumpArrND ExtractReceiver
	cp ./MergeShotAndReceiverASCIIFiles ../../bin
	cp ./DumpArrND ../../bin
	cp ./ExtractReceiver ../../bin

MergeShotAndReceiverASCIIFiles: MergeShotAndReceiverASCIIFiles.cpp
	icc -std=c++11 -I../Geom -I../ArrND -I../../ $^ -L../../lib -lvoxet -o $@

ExtractReceiver: ExtractReceiver.cpp
	icc -std=c++11 -I../Geom -I../ArrND -I../../Common $^ -L../../lib -lvoxet -o $@

DumpArrND: DumpArrND.cpp
	icc -std=c++11 -I../Geom -I../ArrND -I../../Common $^ -L../../lib -lvoxet -o $@
