#include <iostream>
#include <vector>
#include <unordered_map> 
#include <string>
#include <queue>
class Vertex
{
    std::vector<Vertex> neighbours;
    int index, cityIndex;
    double shortestDistance;
    bool visited;
    double x, y;
public:
    Vertex(double x, double y, int cityIndex) : x{x}, y{y}, cityIndex{cityIndex}
    {
        static int vertexAmount;
        visited = false;
        shortestDistance = INFINITY;
        index = vertexAmount++;
    }

    void AddNeighbour(Vertex neighbour)
    {
        neighbours.push_back(neighbour);
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

    double GetShortestDistance()
    {
        return shortestDistance;
    }

    bool operator<(const Vertex& other)
    {
        return this->shortestDistance < other.shortestDistance;
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
    std::cout << "teste";
}

void Djikstras(std::vector<std::string>& cityNames, std::vector<std::vector<Vertex>>& vertexVectors, std::unordered_map<std::string, double>& edgeMap, std::vector<double>& shortestPaths)
{
    std::priority_queue<Vertex> priorityQueue;
    for (int i = 0; i < cityNames.size(); i++)
    {
        shortestPaths.push_back(INFINITY);
    }
    //priorityQueue.push()
    //test priorityqueue and overloaded operator
    //add neihgbors to priorityqueue
    //visit next in queue
    //update distance to vertieses and cities
    //add city-distances together
    //return
}


int main()
{
    std::vector<std::string> cityNames;
    std::vector<std::vector<Vertex>> vertexVectors;
    std::unordered_map<std::string, double> edgeMap;
    std::vector<double> shortestPaths;
    CreateVerticies(cityNames, vertexVectors);
    CreateEdges(cityNames, vertexVectors, edgeMap);
    std::cout << "teste";
}