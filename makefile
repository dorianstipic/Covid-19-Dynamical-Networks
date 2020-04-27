all: model_cluster_trip_v2

model_cluster_trip_v2: model_cluster_trip_v2.cpp
	g++ -O2 -std=c++11 -o model_cluster_trip_v2 model_cluster_trip_v2.cpp
