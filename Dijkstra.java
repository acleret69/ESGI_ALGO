import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.io.IOException;
import java.io.FileReader;
import java.io.BufferedReader;
import java.util.*;

public class Dijkstra {

  // Keep a fast index to nodes in the map
  private Map<String, Vertex> vertexNames;

  /**
   * Construct an empty Dijkstra with a map. The map's key is the name of a vertex
   * and the map's value is the vertex object.
   */
  public Dijkstra() {
    vertexNames = new HashMap<String, Vertex>();
  }

  /**
   * Adds a vertex to the dijkstra. Throws IllegalArgumentException if two vertices
   * with the same name are added.
   * 
   * @param v
   *          (Vertex) vertex to be added to the dijkstra
   */
  public void addVertex(Vertex v) {
    if (vertexNames.containsKey(v.name))
      throw new IllegalArgumentException("Cannot create new vertex with existing name.");
    vertexNames.put(v.name, v);
  }

  /**
   * Gets a collection of all the vertices in the dijkstra
   * 
   * @return (Collection<Vertex>) collection of all the vertices in the dijkstra
   */
  public Collection<Vertex> getVertices() {
    return vertexNames.values();
  }

  /**
   * Gets the vertex object with the given name
   * 
   * @param name
   *          (String) name of the vertex object requested
   * @return (Vertex) vertex object associated with the name
   */
  public Vertex getVertex(String name) {
    return vertexNames.get(name);
  }

  /**
   * Adds a directed edge from vertex u to vertex v
   * 
   * @param nameU
   *          (String) name of vertex u
   * @param nameV
   *          (String) name of vertex v
   * @param cost
   *          (double) cost of the edge between vertex u and v
   */
  public void addEdge(String nameU, String nameV, Double cost) {
    if (!vertexNames.containsKey(nameU))
      throw new IllegalArgumentException(nameU + " does not exist. Cannot create edge.");
    if (!vertexNames.containsKey(nameV))
      throw new IllegalArgumentException(nameV + " does not exist. Cannot create edge.");
    Vertex sourceVertex = vertexNames.get(nameU);
    Vertex targetVertex = vertexNames.get(nameV);
    Edge newEdge = new Edge(sourceVertex, targetVertex, cost);
    sourceVertex.addEdge(newEdge);
  }

  /**
   * Adds an undirected edge between vertex u and vertex v by adding a directed
   * edge from u to v, then a directed edge from v to u
   * 
   * @param nameU
   *          (String) name of vertex u
   * @param nameV
   *          (String) name of vertex v
   * @param cost
   *          (double) cost of the edge between vertex u and v
   */
  public void addUndirectedEdge(String nameU, String nameV, double cost) {
    addEdge(nameU, nameV, cost);
    addEdge(nameV, nameU, cost);
  }


  /**
   * Computes the euclidean distance between two points as described by their
   * coordinates
   * 
   * @param ux
   *          (double) x coordinate of point u
   * @param uy
   *          (double) y coordinate of point u
   * @param vx
   *          (double) x coordinate of point v
   * @param vy
   *          (double) y coordinate of point v
   * @return (double) distance between the two points
   */
  public double computeEuclideanDistance(double ux, double uy, double vx, double vy) {
        // TODO
        double diffx = vx - ux;
        double diffy= vy-uy;
        double sumSquared = Math.pow(diffx, 2) + Math.pow(diffy, 2);
        // return 1.0; // Replace this
        return Math.pow(sumSquared, 0.5);
  }

  /**
   * Calculates the euclidean distance for all edges in the map using the
   * computeEuclideanCost method.
   */
  public void computeAllEuclideanDistances() {
      int dijkstraSize= vertexNames.size();
      for (String u : vertexNames.keySet())
      {
        Vertex node = vertexNames.get(u);
        int sourceX= node.x;
        int sourceY= node.y;
        List adjacencyList=node.adjacentEdges;
        for (Edge e : vertexNames.get(u).adjacentEdges) {

          int targetX=e.target.x;
          int targetY=e.target.y;
          double updatedDistance= computeEuclideanDistance(sourceX, sourceY, targetX, targetY);
          e.distance = updatedDistance;

        }

      }


  }

  
  public void doDijkstra(String s) {
      double MAX_VALUE=1000000;

      for (String u : vertexNames.keySet())
      {
        Vertex node = vertexNames.get(u);
        node.distance=MAX_VALUE;
        node.prev=null;
        node.known=false;

      }
      vertexNames.get(s).distance=0;

      LinkedList<Vertex> vertexList= new LinkedList<>();
          for (String u : vertexNames.keySet())
          {
            vertexList.addLast(vertexNames.get(u));
          }
          
          while(vertexList.size()>0)
          {

            int index= findClosestVertex(vertexList);
            Vertex v= vertexList.get(index);
            v.known=true;
            vertexList.remove(index);

            List<Edge> adjacencyList = v.adjacentEdges;

              for (Edge e : adjacencyList)
              {
                if (!e.target.known)
                {
                  double distanceCost=e.distance;
                  if (v.distance + distanceCost<e.target.distance)
                  {
                    e.target.distance=v.distance + distanceCost;
                    e.target.prev=v;
                  }
                }
              }
          }
        }

           public static int findClosestVertex(LinkedList<Vertex> list)
        {
          double minDist=list.get(0).distance;
          int minIndex=0;
          int i=0;

        for (i=1; i< list.size(); i++)
        { 
          double dist= list.get(i).distance;
          if (dist<minDist)
          {
            minDist=dist;
            minIndex=i;
          }

        }
        return minIndex;
      }

  public List<Edge> getDijkstraPath(String s, String t) {
    doDijkstra(s);
    
    Vertex v = vertexNames.get(t);
    LinkedList<Edge> shortestPathEdges = new LinkedList<>();
    Vertex penultimateVertex = null;

    while(v.prev!=null)
    {
      List<Edge> adjacencyList = v.adjacentEdges; 
      for (Edge e : adjacencyList)
      {
        if (e.target.equals(v.prev))
        {
          shortestPathEdges.addFirst(e);
        }
      }

      penultimateVertex = v;
      v=v.prev;
    }
    for (Edge e : v.adjacentEdges)
    {
        if (e.target.equals(penultimateVertex))
        {
          shortestPathEdges.addFirst(e);
        }
    }    
    return shortestPathEdges; 
  }

  public List<Edge> getdepthCoursePath(String s, String t) {
    LinkedList<Vertex> result = new LinkedList<Vertex>();
    LinkedList<Vertex> summitsVisited = new LinkedList<Vertex>();
    List<Vertex> resultat = doDepthCourse(vertexNames.get(s),result,summitsVisited,vertexNames.get(t));
    LinkedList<Edge> shortestPathEdges = new LinkedList<>();
    for(int i=0; i+1 < resultat.size(); i++) {
      Edge e = new Edge(resultat.get(i), resultat.get(i+1),0);
      shortestPathEdges.add(e);
    }

    return shortestPathEdges;
  }

  public static List<Vertex> doDepthCourse(Vertex T, List<Vertex> result, List<Vertex> summitsVisited,Vertex End) {
    if(summitsVisited.contains(End)){
      return result;
    } else if(!summitsVisited.contains(T)) {
        result.add(T);
        summitsVisited.add(T);
        for(int i=0; i < T.adjacentEdges.size(); i++) {
          doDepthCourse((Vertex) T.adjacentEdges.get(i).target, result, summitsVisited,End);
        }
    }
    return result;
  }
}
