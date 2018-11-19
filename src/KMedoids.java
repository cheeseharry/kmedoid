import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.*;

/**
 * Created by leizheng on 9/29/18.
 */
public class KMedoids {
    public int k;//k cluster

    public static List<List<Double>> origDataList = new ArrayList<>();
    public static List<List<Double>> lastCenterList = new ArrayList<>();
    public static List<List<Double>> centerList = new ArrayList<>();

    //k center with 12 dimension [k][13]
    //public static List<List<Double>> firstCenterList = new ArrayList<>();

    //key=clusetIndex; value=time-series-gene in the cluster
    public static HashMap<Integer,List<List<Double>>> lastMemberShipMap = new HashMap<>();// Store last 10 Cluster
    public static HashMap<Integer,List<List<Double>>> memberShipMap = new HashMap<>();// Store new 10 Cluster

    public static void main(String[] args) {
        KMedoids kmedoids = new KMedoids(10);
        origDataList = getdata();

        //first partition
        centerList = kmedoids.getFirstCenter();

        //get first Cluster
        memberShipMap = kmedoids.calMemberShip2();

        //get new centerList
        lastCenterList = centerList;
        centerList = kmedoids.getNewCenterList(memberShipMap);

        //Loop until no change
        while (!lastMemberShipMap.equals(memberShipMap)){
            //recursive get new membership cluster
            lastMemberShipMap = memberShipMap;
            memberShipMap = kmedoids.calMemberShip2();
            lastCenterList = centerList;
            centerList = kmedoids.getNewCenterList(memberShipMap);
            //System.out.println(centerList);
        }

        //print result
        for (Map.Entry<Integer,List<List<Double>>> entry : memberShipMap.entrySet()){
            List<List<Double>> membershipList = entry.getValue();
            List<Integer> clusterId = new ArrayList<>();
            for (List<Double> geneList : membershipList){
                clusterId.add(geneList.get(0).intValue());
            }
            System.out.println("size "+ membershipList.size() +": " + clusterId);
        }

    }

    /**
     * init
     * @param k
     */
    public KMedoids(int k){
        this.k = k;
    }

    /**
     * get first centerList
     * @return
     */

    public List<List<Double>> getFirstCenter(){
        List<List<Double>> firstCenterList = new ArrayList<>();
        HashMap<List<Double>, Double> firstDistanceMap = getDistanceSumMap(origDataList);
        PriorityQueue<Map.Entry<List<Double>, Double>> firstMinHeap = getMedoidHeap(firstDistanceMap);

        for (int i = 0;i<10;i++){
            Map.Entry entry = firstMinHeap.poll();
            firstCenterList.add((ArrayList)entry.getKey());
        }
        return firstCenterList;
    }

    /**
     * getNewCenterList
     * @param memberShipMap
     * @return
     */
    public List<List<Double>> getNewCenterList(HashMap<Integer,List<List<Double>>> memberShipMap){
        List<List<Double>> newCenterList = new ArrayList<>();
        for (Map.Entry<Integer,List<List<Double>>> entry : memberShipMap.entrySet()){
            int index = entry.getKey();
            List<List<Double>> membershipList = entry.getValue();
            HashMap<List<Double>, Double> nextDistanceMap = getDistanceSumMap(membershipList);
            PriorityQueue<Map.Entry<List<Double>, Double>> nextMinHeap = getMedoidHeap(nextDistanceMap);
            Map.Entry entry2 = nextMinHeap.poll();
            newCenterList.add((ArrayList)entry2.getKey());
        }
        return newCenterList;
    }



    /**
     * Read data from input
     * @return
     */
    public static List<List<Double>> getdata() {
        List<List<Double>> data = new ArrayList<List<Double>>();
        double geneID = 0; //geneID
        try {
            String encoding = "UTF-8";
            File file = new File("assignment3_input.txt");
            if (file.isFile() && file.exists()) {
                InputStreamReader read = new InputStreamReader(
                        new FileInputStream(file), encoding);
                BufferedReader bufferedReader = new BufferedReader(read);
                String lineTXT = null;
                while ((lineTXT = bufferedReader.readLine()) != null) {
                    String[] lineString = lineTXT.split("\\s+");
                    List<Double> lineList = new ArrayList<Double>();
                    for (int i = 0; i < lineString.length; i++) {
                        lineList.add(Double.parseDouble(lineString[i]));
                    }
                    lineList.add(0,geneID);
                    data.add(lineList);
                    geneID++;
                }
                read.close();
            } else {
                System.out.println("file not find！");
            }
        } catch (Exception e) {
            System.out.println("read error");
            e.printStackTrace();
        }
        return data;
    }

    /**
     * Assign cluster for each membership
     */
    public HashMap<Integer,List<List<Double>>> calMemberShip2(){
        HashMap<Integer,List<List<Double>>> tempClusterMap = new HashMap<>();
        for (int j=0; j<origDataList.size();j++){
            double currentDistance = Double.MAX_VALUE;
            int index = -1;
            for (int i=0;i<k;i++){
                List<Double> tempCenter = centerList.get(i);
                //centerList should contain 10 Center with 12 dimension
                double distance = getDistance(origDataList.get(j),centerList.get(i));
                //compare distance with every centerPoint,choose the smallest one
                if (distance<currentDistance){
                    currentDistance = distance;
                    index = i;
                }
            }
            if (tempClusterMap.containsKey(index)){
                tempClusterMap.get(index).add(origDataList.get(j));
            }
            else {
                List<List<Double>> memberShipList = new ArrayList<>();
                memberShipList.add(origDataList.get(j));
                tempClusterMap.put(index,new ArrayList<>(memberShipList));
            }
        }
        return tempClusterMap;
    }

    /**
     * calc Euclidean distance
     * @param gene1
     * @param gene2
     * @return
     */
    public double getDistance(List<Double> gene1, List<Double> gene2){
        double tempDistance = 0;
        if (gene1!=null&&gene2!=null && gene1.size()==gene2.size()){
            //start with index 1 for the 1st dimension
            for (int i=1;i<gene1.size();i++){
                tempDistance += Math.pow(gene1.get(i)-gene2.get(i),2);
            }
        }
        tempDistance = Math.sqrt(tempDistance);
        tempDistance = (double) Math.round(tempDistance * 10000) / 10000;
        return tempDistance;
    }

    /**
     * sumOfDistance
     * @param clusterList
     * @param gene2
     * @return
     */
    public double sumOfDistance(List<List<Double>> clusterList, List<Double> gene2){
        double sumDistance = 0;
        for (List<Double> gene1 : clusterList){
            sumDistance += getDistance(gene1,gene2);
        }
        return sumDistance;
    }

    /**
     * DistanceSum map
     * @param clusterList
     * @return
     */
    public HashMap<List<Double>, Double> getDistanceSumMap(List<List<Double>> clusterList){
        HashMap<List<Double>, Double> distanceMap = new HashMap<>();
        for (List<Double> myList : clusterList){
            distanceMap.put(myList,sumOfDistance(clusterList,myList));
        }

        return distanceMap;
    }

    /**
     * PQ store distanceSum in ascending order
     * @param distanceSumMap
     * @return
     */
    public PriorityQueue<Map.Entry<List<Double>, Double>> getMedoidHeap(HashMap<List<Double>, Double> distanceSumMap){
        PriorityQueue<Map.Entry<List<Double>, Double>> minHeap = new PriorityQueue<>(new Comparator<Map.Entry<List<Double>, Double>>() {
            @Override
            public int compare(Map.Entry<List<Double>, Double> o1, Map.Entry<List<Double>, Double> o2) {
                if (o1.getValue()-o2.getValue() <= 0){
                    return -1;
                }
                else {
                    return 1;
                }
            }
        });
        minHeap.addAll(distanceSumMap.entrySet());

        return minHeap;
    }


}
