all: model_initial model_cluster_trip

model_initial: model_initial.cpp
	g++ -O2 -std=c++11 -o model_initial model_initial.cpp

model_cluster_trip: model_cluster_trip.cpp
	g++ -O2 -std=c++11 -o model_cluster_trip model_cluster_trip.cpp
