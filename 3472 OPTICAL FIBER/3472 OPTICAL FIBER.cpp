#include <iostream>
#include <vector>
#include <unordered_map> 
#include <string>
#include <queue>
class Point
{
public:
    double x, y;

    Point() 
    {
        x = 0;
        y = 0;
    }
};

class Vertex
{
    std::vector<Vertex> neighbours;
    int index, cityIndex;
    double shortestDistance;
    int pathViaIndex;
    bool visited;
    double x, y;
public:
    Vertex(double x, double y, int cityIndex) : x{x}, y{y}, cityIndex{cityIndex}
    {
        static int vertexAmount;
        visited = false;
        shortestDistance = INFINITY;
        pathViaIndex = -1;
        index = vertexAmount++;
    }

    void Visit()
    {
        visited = true;
    }

    bool IsVisited()
    {
        return visited;
    }

    void AddNeighbour(Vertex neighbour)
    {
        neighbours.push_back(neighbour);
    }

    std::vector<Vertex> GetNeighbours()
    {
        return neighbours;
    }

    int GetIndex()
    {
        return index;
    }

    int GetCityIndex()
    {
        return cityIndex;
    }

    double GetX()
    {
        return x;
    }

    double GetY()
    {
        return y;
    }

    void SetPathViaIndex(int index)
    {
        pathViaIndex = index;
    }

    int GetPathViaIndex()
    {
        return pathViaIndex;
    }

    void SetShortestDistance(double distance)
    {
        shortestDistance = distance;
    }

    double GetShortestDistance()
    {
        return shortestDistance;
    }

    friend bool operator<(const Vertex& l, const Vertex& r)
    {
        return l.shortestDistance > r.shortestDistance;
    }
};

class Edge
{
    double distance;
    Vertex start, end;
private:
    double CalculateDistance()
    {
        double x, y, result;
        x = std::abs(end.GetX() - start.GetX());
        y = std::abs(end.GetY() - start.GetY());
        result = sqrt(x * x + y * y);
        return result;
    }
public:
    Edge(Vertex start, Vertex end) : start{start}, end{end}
    {
        distance = CalculateDistance();
    }

    double GetDistance()
    {
        return distance;
    }
};

void CreateVerticies(std::vector<std::string>& cityNames, std::vector<std::vector<Vertex>>& vertexVectors)
{
    int cityAmount = 0;
    int locationAmount = 0;
    int x = 0;
    int y = 0;
    std::string cityName;
    std::vector<Vertex> vertecies;
    std::cin >> cityAmount;
    for (int i = 0; i < cityAmount; i++)
    {
        std::cin >> cityName;
        cityNames.push_back(cityName);
        std::cin >> locationAmount;
        for (int j = 0; j < locationAmount; j++)
        {
            std::cin >> x;
            std::cin >> y;
            Vertex vertex(x, y, i);  //create vertex
            vertecies.push_back(vertex);
        }
        vertexVectors.push_back(vertecies);
        vertecies.clear();
    }
}

std::string GetKey(int startIndex, int endIndex)
{
    return std::to_string(startIndex) + "&" + std::to_string(endIndex);
}

void CreateEdges(std::vector<std::string>& cityNames, std::vector<std::vector<Vertex>>& vertexVectors, std::unordered_map<std::string, double>& edgeMap)
{
    std::string city1, city2;
    int cityIndex1, cityIndex2;
    for (int i = 0; i < cityNames.size() - 1; i++)
    {
        std::cin >> city1;
        std::cin >> city2;
        //Find cities
        for (int j = 0; j < cityNames.size(); j++)
        {
            if (cityNames.at(j) == city1)
            {
                cityIndex1 = j;
                break;
            }
        }
        for (int j = 0; j < cityNames.size(); j++)
        {
            if (cityNames.at(j) == city2)
            {
                cityIndex2 = j;
                break;
            }
        }
        
        //Create edges
        for (Vertex v : vertexVectors.at(cityIndex1))
        {
            for (Vertex w : vertexVectors.at(cityIndex2))
            {
                v.AddNeighbour(w);
                w.AddNeighbour(v);
                Edge edge(v, w);
                edgeMap[GetKey(v.GetIndex(), w.GetIndex())] = edge.GetDistance();
                edgeMap[GetKey(w.GetIndex(), v.GetIndex())] = edge.GetDistance();
            }
        }
    }
}

//void Djikstras(std::vector<std::string>& cityNames, std::vector<std::vector<Vertex>>& vertexVectors, std::unordered_map<std::string, double>& edgeMap, std::vector<double>& shortestCityPaths)
//{
//    std::priority_queue<Vertex&> priorityQueue;
//    for (int i = 0; i < cityNames.size(); i++)
//    {
//        shortestCityPaths.push_back(INFINITY);
//    }
//
//    //Create starting vertex
//    Vertex start(0, 0, 0);
//    for (int i = 0; i < vertexVectors.at(0).size(); i++)
//    {
//        edgeMap[GetKey(start.GetIndex(), vertexVectors.at(0).at(i).GetIndex())] = 0;
//        edgeMap[GetKey(vertexVectors.at(0).at(i).GetIndex(), start.GetIndex())] = 0;
//    }
//    start.SetShortestDistance(0);
//    priorityQueue.push(start);
//
//    while (!priorityQueue.empty())
//    {
//        Vertex& current = priorityQueue.top();
//        priorityQueue.pop();
//        current.Visit();
//        for (Vertex v : current.GetNeighbours()) 
//        {
//            if (v.GetShortestDistance() > current.GetShortestDistance() + edgeMap[GetKey(current.GetIndex(), v.GetIndex())])
//            {
//                v.SetShortestDistance(current.GetShortestDistance() + edgeMap[GetKey(current.GetIndex(), v.GetIndex())]);
//                v.SetPathViaIndex(current.GetIndex());
//                if (v.GetShortestDistance() < shortestCityPaths.at(v.GetCityIndex())) shortestCityPaths.at(v.GetCityIndex()) = v.GetShortestDistance();
//            }
//        }
//    }
//
//
//    //add neihgbors to priorityqueue
//    //visit next in queue
//    //update distance to vertieses and cities
//    //add city-distances together
//    //return
//}

Point GetCityAverage(std::vector<std::string>& cityNames, std::vector<std::vector<Vertex>>& vertexVectors)
{
    Point cityAverage;
    std::vector<Point> averageList;
    for (int i = 0; i < cityNames.size(); i++)
    {
        Point average;
        for (int j = 0; j < vertexVectors.at(i).size(); j++)
        {
            average.x += vertexVectors.at(i).at(j).GetX();
            average.y += vertexVectors.at(i).at(j).GetY();
        }
        average.x = average.x / vertexVectors.at(i).size();
        average.y = average.y / vertexVectors.at(i).size();
        averageList.push_back(average);
    }
    for (int i = 0; i < averageList.size(); i++)
    {
        cityAverage.x += averageList.at(i).x;
        cityAverage.y += averageList.at(i).y;
    }
    cityAverage.x = cityAverage.x / averageList.size();
    cityAverage.y = cityAverage.y / averageList.size();

    return cityAverage;
}

double ShortestPath(std::vector<std::string>& cityNames, std::vector<std::vector<Vertex>>& vertexVectors, std::unordered_map<std::string, double>& edgeMap)
{
    double shortestPath = INFINITY;
    Point average = GetCityAverage(cityNames, vertexVectors);
    return shortestPath;
}


int main()
{
    std::vector<std::string> cityNames;
    std::vector<std::vector<Vertex>> vertexVectors;
    std::unordered_map<std::string, double> edgeMap;
    std::vector<double> shortestCityPaths;
    CreateVerticies(cityNames, vertexVectors);
    CreateEdges(cityNames, vertexVectors, edgeMap);
    ShortestPath(cityNames, vertexVectors, edgeMap);
}