using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Text;
using System.Windows.Forms;
using System.IO;
using System.Media;
namespace ImageQuantization
{
    public partial class MainForm : Form
    {
        public MainForm()
        {
            InitializeComponent();
        }
        double totalweight;
        RGBPixel[,] ImageMatrix;

        private void btnOpen_Click(object sender, EventArgs e)
        {
            OpenFileDialog openFileDialog1 = new OpenFileDialog();
            if (openFileDialog1.ShowDialog() == DialogResult.OK)
            {
                //Open the browsed image and display it
                string OpenedFilePath = openFileDialog1.FileName;
                ImageMatrix = ImageOperations.OpenImage(OpenedFilePath);
                ImageOperations.DisplayImage(ImageMatrix, pictureBox1);
            }
            txtWidth.Text = ImageOperations.GetWidth(ImageMatrix).ToString();
            txtHeight.Text = ImageOperations.GetHeight(ImageMatrix).ToString();
        }

        private void btnGaussSmooth_Click(object sender, EventArgs e)
        {
            double sigma = double.Parse(txtGaussSigma.Text);
            int maskSize = (int)nudMaskSize.Value ;
            ImageMatrix = ImageOperations.GaussianFilter1D(ImageMatrix, maskSize, sigma);
            ImageOperations.DisplayImage(ImageMatrix, pictureBox2);
        }

        private void button1_Click(object sender, EventArgs e)
        {
            
            string temp = "";
            Node mynode = new Node();
            Graph G = new Graph();
            Dictionary<string, int> mydic = new Dictionary<string, int>();
            int mi = 0;
          
            for (int i = 0; i < ImageOperations.GetHeight(ImageMatrix); i++)
            {
                for (int j = 0; j < ImageOperations.GetWidth(ImageMatrix); j++)
                {
                    temp += ImageMatrix[i, j].red.ToString() + ','
                            + ImageMatrix[i, j].green.ToString() + ','
                            + ImageMatrix[i, j].blue.ToString();
                    if (!mydic.ContainsKey(temp))
                    {
                        
                        mydic.Add(temp, mi);
                        mynode = new Node();
                        mynode.name = mi;
                        mynode.val = ImageMatrix[i, j];
                        G.AddNodeToGraph(mynode);
                        mi++;
                    }
                    temp = "";
                }
            }

            int k = 0;
            List<KeyValuePair<double, KeyValuePair<int, int>>> nt=null;
            if (string.IsNullOrEmpty(textBox1.Text))
            {
                List<KeyValuePair<double, KeyValuePair<int, int>>> myTree = MST(G);
                List<KeyValuePair<double, KeyValuePair<int, int>>> newT = new List<KeyValuePair<double, KeyValuePair<int, int>>>(myTree);
                nt = newT;
                 k = clusternum(myTree, G);
                 textBox1.Text= k.ToString();
            }
            else {
                k = int.Parse(textBox1.Text);
                
            }
            List<RGBPixel> centroids;
            if (comboBox1.SelectedItem.ToString() == "SingleLinkage")
            {
                if(nt==null)
                    nt = MST(G);
              List<set> mysets=  newclustering(nt, k, G);
             centroids = setcolors(mysets, ref G);
            }
            else if (comboBox1.SelectedItem.ToString() == "HEMST")
            {
                List<set> a = clusterbonus(G, k);
                 centroids = setcolors(a, ref G);
            }
            else {
                Console.WriteLine(G.GraphList.Count);
                splittreeNode tha = splittree(G.GraphList);
                List<KeyValuePair<splittreeNode, splittreeNode>> p = towspd(tha);
                nt = makemst(p, G);
                List<set> returned = wspdclustering(nt, k);
                centroids = setcolors(returned, ref G);
            }

                temp = "";
            Console.WriteLine(centroids.Count);
                RGBPixel[,] narr = new RGBPixel[ImageOperations.GetHeight(ImageMatrix), ImageOperations.GetWidth(ImageMatrix)];
                for (int i = 0; i < ImageOperations.GetHeight(ImageMatrix); i++)
                {
                    for (int j = 0; j < ImageOperations.GetWidth(ImageMatrix); j++)
                    {
                        temp += ImageMatrix[i, j].red.ToString() + ','
                                + ImageMatrix[i, j].green.ToString() + ','
                                + ImageMatrix[i, j].blue.ToString();

                        narr[i, j] = centroids[G.GraphList[mydic[temp]].centnum];                 
                        temp = "";
                    }
                }
                ImageOperations.DisplayImage(narr, pictureBox2);
            }
       KeyValuePair<double, KeyValuePair<set, int>> getEdge(double bigedge,set s,Graph G) {
            Queue<set> myque = new Queue<set>();
            set parent= new set();
            int index = -1;
            double maxdiff = -1;
            myque.Enqueue(s);
            while (myque.Count != 0) {
                for (int i = 0; i < myque.Peek().children.Count; i++) {
                   double mymin=Math.Min(s.count - myque.Peek().children[i].count, myque.Peek().children[i].count);
                    if (maxdiff < mymin &&
                        G.getcost(myque.Peek().children[i].parent.name, myque.Peek().children[i].name) == bigedge)
                    {

                            maxdiff = mymin;
                            parent = myque.Peek();
                            index = i;
                    
                    }


                    myque.Enqueue(myque.Peek().children[i]);
                }
                myque.Dequeue();
            
            }
            KeyValuePair<set, int> myvals = new KeyValuePair<set, int>(parent,index);
            KeyValuePair<double, KeyValuePair<set, int>> k = new KeyValuePair<double, KeyValuePair<set, int>>(maxdiff,myvals);
            return k;
        
        
        }
        List<set> newclustering(List<KeyValuePair<double, KeyValuePair<int, int>>> T,int k,Graph G) {
            T.Sort((x, y) => x.Key.CompareTo(y.Key));  // uses quick sorting complexity is(nlogn)
            List<set> myset= makedisjoint(T);          // returns a set of each node in the graph combined in a disjoint set
            List<set> clusters = new List<set>();      // complexity (it traverses to all the nodes) n
            int l = T.Count - 1;                      // gets the parent of the whole cluster normal time takes (log n) and n
            set parent = findparent(myset[0]);         // in the worst case
            treecounting(ref parent,G);                //sets standard deviation mean and number of nodes using bfs 
            clusters.Add(parent);                      //every node will be visited therefore n

            while (clusters.Count != k) {             // each statement will execute k times
                double value = T[l].Key;
                int index = -1;
                set todelete = new set();
               set myparent= new set();
                l--;
                double TotalMax=0;
                for (int i = 0; i < clusters.Count; i++)   //starts with 1 iteration the increases in value by 1
                {                                         //1+2+3+4+......+k-1=k*(k-1)/2 therefore (k^2) every statement except
                  KeyValuePair<double,KeyValuePair<set, int>> pair = getEdge(value, clusters[i], G);  // getedge
                  double max = pair.Key;                            // it uses bfs algorithm to traverse through
                  if (max > TotalMax) {                            // all the nodes of a given cluster therefore total
                      todelete = pair.Value.Key;                  // time it takes evey iteration is n and the above iteration
                       TotalMax=max;                             // complexity is k therefore kn
                      index = pair.Value.Value;
                      myparent = clusters[i];
                  
                  }
                }
                todelete.children[index].parent = null;
                set child = todelete.children[index];
                todelete.children[index] = todelete.children[todelete.children.Count - 1];
                todelete.children.RemoveAt(todelete.children.Count - 1);
                treecounting(ref parent,G);
                clusters.Add(child);
            }

            return clusters;
        }


        int treecounting(ref set s,Graph G) {
            s.Edgesvalue = 0;
            s.sqrEdgesvalue = 0;
            if (s.children.Count == 0)
            {
                s.count = 1;
                return 1;
               
            }
            int childcount = 0;
            double Edgevals = 0;
            double sqrEdgevals = 0;
            for (int i = 0; i < s.children.Count; i++) {
                set myset = s.children[i];
                childcount += treecounting(ref myset,G);
                Edgevals += G.getcost(s.children[i].name, s.name) + s.children[i].Edgesvalue;
                sqrEdgevals += G.getcost(s.children[i].name, s.name) * G.getcost(s.children[i].name, s.name)+s.children[i].sqrEdgesvalue;

            }
            s.Edgesvalue = Edgevals;
            s.sqrEdgesvalue = sqrEdgevals;
            s.count = 1 + childcount;
            return s.count;
        
        
        }

        List<KeyValuePair<double, KeyValuePair<int, int>>> makemst(List<KeyValuePair<splittreeNode,splittreeNode>> p,Graph G) {
            List<KeyValuePair<double, KeyValuePair<int, int>>> mylist = new List<KeyValuePair<double, KeyValuePair<int, int>>>();
            double mynum = 0;
            for (int i = 0; i < p.Count; i++) {
                double cost = G.getcost(p[i].Key.num, p[i].Value.num);
                KeyValuePair<int,int> k= new KeyValuePair<int,int>(p[i].Key.num, p[i].Value.num);
                KeyValuePair<double, KeyValuePair<int, int>> k2 = new KeyValuePair<double, KeyValuePair<int, int>>(cost,k);
                mylist.Add( k2);
                mynum += cost;
            }
            return mylist;
        }
      
        
        splittreeNode splittree(List<Node> G)
        {
            splittreeNode mytree = new splittreeNode();
            if (G.Count == 1)
            {
                splittreeNode mynode = new splittreeNode();
                mynode.val = G[0].val;
                mynode.num = G[0].name;
                mynode.avg = G[0].val;
                return mynode;
            }
            byte mingreen = 255;
            byte minred = 255;
            byte minblue = 255;
            byte maxgreen = 0;
            byte maxred = 0;
            byte maxblue = 0;
            for (int i = 0; i < G.Count; i++)        //gets bounds of the cube n complexity 
            {

                if (G[i].val.green < mingreen)
                {
                    mingreen = G[i].val.green;
                }
                if (G[i].val.red < minred)
                {
                    minred = G[i].val.red;
                }
                if (G[i].val.blue < minblue)
                {
                    minblue = G[i].val.blue;
                }

                if (G[i].val.green > maxgreen)
                {
                    maxgreen = G[i].val.green;
                }
                if (G[i].val.red > maxred)
                {
                    maxred = G[i].val.red;
                }
                if (G[i].val.blue > maxblue)
                {
                    maxblue = G[i].val.blue;
                }
            }
            byte diffred = (byte)(maxred - minred);
            byte diffgreen = (byte)(maxgreen - mingreen);
            byte diffblue = (byte)(maxblue - minblue);
            int choice = 0;
            if (diffred >= diffgreen && diffred >= diffblue)
            {
                choice = 0;
            }
            else if (diffgreen >= diffred && diffgreen >= diffblue)
            {
                choice = 1;
            }
            else
                choice = 2;
            List<Node> leftnodes = new List<Node>();
            List<Node> rightnodes = new List<Node>();
            byte redavg = (byte)((maxred + minred) / 2);
            byte greenavg = (byte)((maxgreen + mingreen) / 2);
            byte blueavg = (byte)((maxblue + minblue) / 2);

            if (choice == 0)
            {
                byte avg = (byte)((maxred + minred) / 2);
                for (int i = 0; i < G.Count; i++)
                {
                    if (G[i].val.red <= avg)
                    {
                        leftnodes.Add(G[i]);
                    }
                    else
                        rightnodes.Add(G[i]);
                }

            }
            if (choice == 1)
            {
                byte avg = (byte)((maxgreen + mingreen) / 2);
                for (int i = 0; i < G.Count; i++)
                {
                    if (G[i].val.green <= avg)
                    {
                        leftnodes.Add(G[i]);
                    }
                    else
                        rightnodes.Add(G[i]);
                }

            }
            if (choice == 2)
            {
                byte avg = (byte)((maxblue + minblue) / 2);
                for (int i = 0; i < G.Count; i++)
                {
                    if (G[i].val.blue <= avg)
                    {
                        leftnodes.Add(G[i]);
                    }
                    else
                        rightnodes.Add(G[i]);
                }

            }
           
            splittreeNode myNode = new splittreeNode();
            myNode.greenvals = new KeyValuePair<byte, byte>(mingreen, maxgreen);
            myNode.bluevals = new KeyValuePair<byte, byte>(minblue, maxblue);
            myNode.redvals = new KeyValuePair<byte, byte>(minred, maxred);
            myNode.avg.red = redavg;
            myNode.avg.green = greenavg;
            myNode.avg.blue = blueavg;
            myNode.left = splittree(leftnodes);
            myNode.right = splittree(rightnodes);
            return myNode;

        }
        List<KeyValuePair<splittreeNode,splittreeNode>>  towspd(splittreeNode T) {
            List<KeyValuePair<splittreeNode, splittreeNode>> mylist = new List<KeyValuePair<splittreeNode, splittreeNode>>();
            Queue<splittreeNode> nodes = new Queue<splittreeNode>();
            nodes.Enqueue(T);
            while (nodes.Count != 0) {                        //for each non leaf node in splittree it generates a unique pair 
                if (nodes.Peek().num == -1)
                {
                    nodes.Enqueue(nodes.Peek().left);
                    nodes.Enqueue(nodes.Peek().right);
                    mylist.Add(findpair(nodes.Peek().left, nodes.Peek().right));
                }
                nodes.Dequeue();
            }

            return mylist;
        }

        KeyValuePair<splittreeNode, splittreeNode> findpair(splittreeNode u,splittreeNode v) {
            KeyValuePair<splittreeNode, splittreeNode> mypair;
            splittreeNode first = new splittreeNode();
            splittreeNode second = new splittreeNode();
            if (u.num != -1 && v.num != -1) { 
            mypair= new KeyValuePair<splittreeNode,splittreeNode>(u,v);
                return mypair;
            }
            else if (u.num == -1) { 
            first = findpair(u.left, v).Key;
            second = findpair(u.right, v).Value;
            mypair = new KeyValuePair<splittreeNode, splittreeNode>(first, second);
            return mypair;
            }
            else if (v.num == -1) {
                first = findpair(u, v.left).Key;
                second = findpair(u, v.right).Value;
                mypair = new KeyValuePair<splittreeNode, splittreeNode>(first, second);
                return mypair;
            }
            if (getmax(u) < getmax(v))
            {
                first = findpair(u.left, v).Key;
                second = findpair(u.right, v).Value;
            }
            else {
                first = findpair(u, v.left).Key;
                second = findpair(u, v.right).Value;
            }
        mypair = new KeyValuePair<splittreeNode, splittreeNode>(first,second);
        return mypair;
        }
        int getmax(splittreeNode u) {
            return Math.Max(Math.Max(u.greenvals.Value - u.greenvals.Key, u.redvals.Value - u.redvals.Key), u.bluevals.Value - u.bluevals.Key);
        }
        int clusternum(List<KeyValuePair<double, KeyValuePair<int, int>>> Tree, Graph G) {
            int k = 0;
            List<set> s = makedisjoint(Tree);       //n
            set parent = findparent(s[0]);          //logn
            List<set> clusters = new List<set>();
            treecounting(ref parent, G);           //n
            double std = parent.sqrEdgesvalue/(parent.count-1);
            int count = parent.count - 1;
            std -= (parent.Edgesvalue / (parent.count - 1))*(parent.Edgesvalue / (parent.count - 1));
            if (std < 0.001 && std > -0.001)
                std = 0;
            std = Math.Sqrt(std);
            if (std == 0)
                return 1;
            clusters.Add(parent);
            List<double> stdlist = new List<double>();
            KeyValuePair<double, double> mypair = new KeyValuePair<double, double>();
            while(true){
                double tempstd = 0;
                set myparent;
                int first = -1;
                int second = -1;
                int finalfirst = -1;
                int finalsecond = -1;
                int deleteindex = 0;
                double nxt = 0; 
                for (int i = 0; i<Tree.Count; i++)         //n
              {
                   first = Tree[i].Value.Key;
                   second = Tree[i].Value.Value;
                   myparent = findparent(s[first]);                 //logn
                   clusters.Add(s[second]);
                   myparent.Edgesvalue -= s[second].Edgesvalue;
                   myparent.Edgesvalue -= Tree[i].Key;
                   myparent.sqrEdgesvalue -= s[second].sqrEdgesvalue ;
                   myparent.sqrEdgesvalue -= Tree[i].Key*Tree[i].Key;
                   myparent.count -= (s[second].count);
                   double avgstd = 0;
                   for (int j = 0; j < clusters.Count; j++) {      // k therefore kn in each outer loop itertion
                       if (clusters[j].count != 1 )
                       {
                           double val = clusters[j].sqrEdgesvalue / (clusters[j].count - 1);
                           val -= (clusters[j].Edgesvalue / (clusters[j].count - 1)) * (clusters[j].Edgesvalue / (clusters[j].count - 1));
                           if (val < 0.001 && val > -0.001)
                               val = 0;
                           else
                           {
                               val = Math.Sqrt(val);
                               val = val * (clusters[j].count - 1);
                           }
                           avgstd += val;
                       }
                   }
                   avgstd /= count-1;
                   
                    if (tempstd < std - avgstd)
                   {
                       tempstd = Math.Abs(avgstd-std);
                       finalfirst = first;
                       finalsecond = second;
                       deleteindex = i;
                       nxt = avgstd;
                   }
                   clusters.RemoveAt(clusters.Count-1);
                   myparent.Edgesvalue += s[second].Edgesvalue + Tree[i].Key;
                   myparent.sqrEdgesvalue += s[second].sqrEdgesvalue + Tree[i].Key*Tree[i].Key;
                   myparent.count += s[second].count;
                
               }
                    int index = s[finalfirst].children.IndexOf(s[finalsecond]); //parent has small children can be considered constant 
                    s[finalfirst].children.RemoveAt(index);                     // by mst structure
                    s[finalsecond].parent = null;
                    clusters.Add(s[finalsecond]);
                    stdlist.Add(tempstd);
                
                if(stdlist.Count>2)
                    if (Math.Abs(stdlist[k] - stdlist[k - 1]) < 0.0001)
                {
                       return k+1;
                }
                k++;
                
                myparent = findparent(s[finalfirst]);              //log n
                treecounting(ref myparent,G);                      //n
                Tree[deleteindex] = Tree[Tree.Count - 1];
                Tree.RemoveAt(Tree.Count - 1);
                if (nxt == 0)
                    return k+1;
                std =nxt;
                count--;
            }
        }
        public List<KeyValuePair<double, KeyValuePair<int, int>>> MST(Graph G)
        {         // minimum spanning tree using prim algorithm
           List< KeyValuePair<double, KeyValuePair<int, int>> >totaltree = new List<KeyValuePair<double,KeyValuePair<int,int>>> (G.GraphList.Count);
            regulator graphtree = new regulator(G.GraphList.Count);
            double mynum=0;
            KeyValuePair<int, int> pair;
            double min=500;
            int Nodenum = -1;
            int other = -1;
            int second=0;
            int index=-1;                                    // it loops through v therefore theta(v)
            for (int i = 1; i < G.GraphList.Count; i++) {   //loops through all verticies except 0 and calculates all costs of
            double cost = G.getcost(0,i);                   // nodes connected to zero
            graphtree.insert(cost, 0, i);                   
                if(cost<min){                               // since we are traversing through all the graph min edge can be
                    min = cost;                             // calculated during looping
                    other = 0;
                    Nodenum = i;
                    index = i-1;
                }
            }
            mynum += min;

            graphtree.deletemin(index);    // the edge is deleted with the node and inserted into the mst
            pair = new KeyValuePair<int, int>(other, Nodenum);
            KeyValuePair<double, KeyValuePair<int, int>> thepair = new KeyValuePair<double, KeyValuePair<int, int>>(min, pair);
            totaltree.Add(thepair);
            
            while(!graphtree.empty())                       //number of graph nodes - 1 loops v-1 times
            {
                min = 500;

                for (int j = 0; 
                    j < graphtree.myList.Count; j++)       // count = v-1 in the first time and each time it dec by 1
                {                                          // therefore v-1+v-2+v-3... sequence is the complexity
                                                           //  (v-1)*(v)/2 = E
                    double cost = G.getcost(Nodenum, graphtree.pairvals[j].Value);
                   
                    if (graphtree.myList[j] > cost)         // each time calculates all the cost with the remaining items and checks
                    {                                       //if it's the minimum node connected to each item
                        graphtree.myList[j] = cost;         
                        pair = new KeyValuePair<int, int>(Nodenum, graphtree.pairvals[j].Value);
                        graphtree.pairvals[j] = pair;        
                    }
                    else
                        cost = graphtree.myList[j];
                   
                   if (cost < min)                             //traversing through the graphgets us the min edge
                     {
                            min = cost;
                            other = graphtree.pairvals[j].Key;
                            second = graphtree.pairvals[j].Value;
                            index = j;
                     }
                   
                }
                pair = new KeyValuePair<int, int>(other, second);
                thepair = new KeyValuePair<double, KeyValuePair<int, int>>(min, pair);
                totaltree.Add(thepair);
                graphtree.deletemin(index);
                Nodenum = second;
                mynum += min;
              
            }
           // MessageBox.Show(mynum.ToString());
            return totaltree;
        }
      List<set> wspdclustering(List<KeyValuePair<double,KeyValuePair<int,int>>> T,int k) {
            List<set> myset = new List<set>(T.Count);
            List<set> clusters = new List<set>();
            for (int i = 0;i< T.Count+1; i++) {           // O(n)
                
                set aset = new set();
                aset.name = i;
                myset.Add(aset);

                clusters.Add(myset[i]);
                clusters[i].clusterindex = i;
            }
            int index = 0;
            T.Sort((x, y) => x.Key.CompareTo(y.Key));     //nlogn
            int l = 0;
            while (k!=clusters.Count) {                   //worst case k=1 will execute n times
                int first = T[l].Value.Key;
                int second = T[l].Value.Value;
                l++;
               set fset=findparent(myset[first]);        //log n 
               set sset = findparent(myset[second]);

               if (fset.root > sset.root)
                   {
                       myset[fset.name].children.Add(myset[sset.name]);
                       myset[sset.name].parent = myset[fset.name];
                       
                        index = myset[sset.name].clusterindex;
                        clusters[myset[fset.name].clusterindex] = myset[fset.name];
                      if (index != clusters.Count - 1)
                       {
                       clusters[clusters.Count - 1].clusterindex = index;
                        clusters[index] = null;
                        clusters[index] = clusters[clusters.Count - 1];
                      }
                      clusters.RemoveAt(clusters.Count - 1);
                   }
               else if (fset.root < sset.root)
                   {
                       myset[sset.name].children.Add(myset[fset.name]);
                       myset[fset.name].parent = myset[sset.name];
                       
                       index = myset[fset.name].clusterindex;
                       clusters[myset[sset.name].clusterindex] = myset[sset.name];
                       if (index != clusters.Count - 1)
                       {
                           clusters[index] = null;
                           clusters[index] = clusters[clusters.Count - 1];

                           clusters[clusters.Count - 1].clusterindex = index;
                       }
                       clusters.RemoveAt(clusters.Count - 1);
                       }
                   else
                   {
                       myset[fset.name].root++;
                       myset[fset.name].children.Add(myset[sset.name]);
                       myset[sset.name].parent = myset[fset.name];
                       index = myset[sset.name].clusterindex;
                       clusters[myset[fset.name].clusterindex] = myset[fset.name];
                       if (index != clusters.Count - 1)
                       {
                           clusters[index] = null;
                           clusters[index] = clusters[clusters.Count - 1];
                           clusters[index].clusterindex = index;
                       }
                       clusters.RemoveAt(clusters.Count - 1);
               }
              
            }
            
           
            return clusters;
        }
        public set findparent(set s)
        {
            set myset = new set();
            if (s.parent != null)
            {
              return findparent(s.parent);

            }
            myset = s;
            return myset;
        }

        public List<RGBPixel> setcolors(List<set> L, ref Graph G)
        {
            List<RGBPixel> centroids = new List<RGBPixel>(L.Count);
            for (int i = 0; i < L.Count; i++) {              // in each cluster it calculates the centroids
                centroids.Add(traverse(L[i], ref G, i));
            }

                return centroids;
        }
        public RGBPixel traverse(set s, ref Graph G,int k) {
            int count=0;
            long bluevals = 0;
            long redvals = 0;
            long greenvals = 0;
            Queue<set> myqueue = new Queue<set>();
            myqueue.Enqueue(s);
            while (myqueue.Count != 0) {                                  //using bfs algorithm each cluster is traversed in
                for (int i = 0; i < myqueue.Peek().children.Count; i++) { //complexity d(d is the number of nodes in this cluster)
                    myqueue.Enqueue(myqueue.Peek().children[i]);          //therefore for the whole algorithim it takes O(n)
                }                                                        //as it's the summation of nodes in each cluster
                int node = myqueue.Peek().name;
                G.GraphList[node].centnum = k;
                bluevals += G.GraphList[node].val.blue;
                redvals += G.GraphList[node].val.red;
                greenvals += G.GraphList[node].val.green;
                count++;
                myqueue.Dequeue();
            }
            bluevals /= count;
            redvals /= count;
            greenvals /= count;
            RGBPixel R = new RGBPixel();
            R.red = (byte) redvals;
            R.blue = (byte)bluevals;
            R.green = (byte)greenvals;
            return R;
     
        }
       public List<set> clusterbonus(Graph G,int k) {
           List<set> clusters = new List<set>();
           List<set> prevclusters = null;
           while (true)
           {
              List<KeyValuePair<double, KeyValuePair<int, int>>> mylist = new List<KeyValuePair<double, KeyValuePair<int, int>>>();
              mylist = MST(G);                                   //E
              mylist.Sort((x, y) => x.Key.CompareTo(y.Key));     //nlogn
              List<set> myset = makedisjoint(mylist);            //n
              double mean = 0;
              for (int i=0; i<mylist.Count;i++ ) {               //n
                  mean += mylist[i].Key;
              }
             int  count=mylist.Count;
             double std = 0;
              mean /= (count);
              for (int i = 0; i < mylist.Count; i++) {          //n
                  std += mylist[i].Key * mylist[i].Key;
              }
              std /= count;
              std -= mean *mean;
              std = Math.Sqrt(std);
              set bset = myset[0];
              treecounting(ref bset,G);                         //n
              bset.clusterindex = 0;
              clusters.Add(bset);
              
              while (mylist[mylist.Count - 1].Key>std+mean) {   //edge cut
                  int num =mylist[mylist.Count - 1].Value.Key;
                  int num2 = mylist[mylist.Count - 1].Value.Value;
                  myset[num].children.RemoveAt(myset[num].children.Count - 1);
                  myset[num2].parent = null;
                  clusters.Add(myset[num2]);
                  mylist.RemoveAt(mylist.Count - 1);
              }
              for (int i = 0; i < clusters.Count; i++) {         //n
                  set theset = clusters[i];
                  treecounting(ref theset, G);
              
              }
                  if (clusters.Count <= k)                     
                  {
                      while (clusters.Count < k)
                      {
                          double value = mylist[mylist.Count - 1].Key;
                          int index = -1;
                          set todelete = new set();
                          set myparent = new set();

                          double TotalMax = 0;
                          for (int i = 0; i < clusters.Count; i++)
                          {
                              KeyValuePair<double, KeyValuePair<set, int>> pair = getEdge(value, clusters[i], G);
                              double max = pair.Key;
                              if (max > TotalMax)
                              {
                                  todelete = pair.Value.Key;
                                  TotalMax = max;
                                  index = pair.Value.Value;
                                  myparent = clusters[i];
                              }
                          }
                          todelete.children[index].parent = null;
                          set child = todelete.children[index];
                          todelete.children[index] = todelete.children[todelete.children.Count - 1];
                          todelete.children.RemoveAt(todelete.children.Count - 1);
                          treecounting(ref myparent, G);
                          clusters.Add(child);
                          mylist.RemoveAt(mylist.Count - 1);

                      }
                      if (prevclusters != null)
                      {
                          List<set> oclusters = new List<set>();
                          for (int i = 0; i < clusters.Count; i++)
                          {
                              oclusters.Add(merger(clusters[i], prevclusters, i));
                          }
                          for (int i = 0; i < oclusters.Count; i++)
                          {
                              oclusters[i].clusterindex = i;
                          }
                          prevclusters = oclusters;

                          return prevclusters;
                      }
                      else
                      {
                          return clusters;
                      }
                  }
                  else
                  {
                      if (prevclusters != null)
                      {
                          List<set> oclusters = new List<set>();
                          for (int i = 0; i < clusters.Count; i++)              //merges two clusters in linear time
                          {
                              oclusters.Add(merger(clusters[i], prevclusters, i));
                          }
                          for (int i = 0; i < oclusters.Count; i++)            // give clusters thier names in linear time
                          {
                              oclusters[i].clusterindex = i;
                          }
                          prevclusters = oclusters;
                      }
                      Graph myG = new Graph();
                      List<set> aset = new List<set>();
                      for (int i = 0; i < clusters.Count; i++)    //traverse each node in each cluster to get avg and min
                      {                                           //using bfs search complexity n
                          Node n = new Node();                    
                          n = traversemin(clusters[i], G, traverseavg(clusters[i], G));
                          myG.AddNodeToGraph(n);
                          set s = new set();
                          s.name = n.name;
                          aset.Add(s);
                      }
                      if (prevclusters == null)
                          prevclusters = clusters;
                      myset = aset;
                      G = myG;
                      clusters = new List<set>();
                  }
           }
        }
        public set merger(set nclusters,List<set> oclusters,int k){
            set s=new set();
            set res = new set();
            Queue<set> myqueue = new Queue<set>();
            myqueue.Enqueue(nclusters);
            int n = 0;
            while (myqueue.Count != 0)
            {
                for (int i = 0; i < myqueue.Peek().children.Count; i++)
                {
                    myqueue.Enqueue(myqueue.Peek().children[i]);
                }
                res = myqueue.Peek();
                int node = myqueue.Peek().name;
                if (n == 0)
                {
                    s = oclusters[node];
                    s.clusterindex = k;
                    n++;
                }
                else {
                    oclusters[node].clusterindex = -1;
                    s.children.Add(oclusters[node]);
                    
                }
                myqueue.Dequeue();
            }

            return s;
        }

       public RGBPixel traverseavg(set s, Graph G){
           int count = 0;
           long bluevals = 0;
           long redvals = 0;
           long greenvals = 0;
           Queue<set> myqueue = new Queue<set>();
           myqueue.Enqueue(s);
           while (myqueue.Count != 0)
           {
               for (int i = 0; i < myqueue.Peek().children.Count; i++)
               {
                   myqueue.Enqueue(myqueue.Peek().children[i]);
               }
               int node = myqueue.Peek().name;
               bluevals += G.GraphList[node].val.blue;
               redvals += G.GraphList[node].val.red;
               greenvals += G.GraphList[node].val.green;
               count++;
               myqueue.Dequeue();
           }
           bluevals /= count;
           redvals /= count;
           greenvals /= count;
           RGBPixel R = new RGBPixel();
           R.red = (byte)redvals;
           R.blue = (byte)bluevals;
           R.green = (byte)greenvals;
           return R;

       }
       public Node traversemin(set s, Graph G,RGBPixel avg)
       {
           Node pic = new Node();
           Queue<set> myqueue = new Queue<set>();
           myqueue.Enqueue(s);
           double min = 500;
           while (myqueue.Count != 0)
           {
               for (int i = 0; i < myqueue.Peek().children.Count; i++)
               {
                   myqueue.Enqueue(myqueue.Peek().children[i]);
               }
               int node = myqueue.Peek().name;
               double red = G.GraphList[node].val.red - avg.red;
               red *= red;
               double green = G.GraphList[node].val.green - avg.green;
               green *= green;
               double blue = G.GraphList[node].val.blue - avg.blue;
               blue *= blue;
               double cost = Math.Sqrt(red + green + blue);
               if (cost < min) {
                   min = cost;
                   pic = G.GraphList[node];
               }
               myqueue.Dequeue();
           }
           return pic;

       }
      List<set> makedisjoint(List<KeyValuePair<double,KeyValuePair<int,int>>> T) {
           List<set> myset = new List<set>();
           for (int i = 0; i < T.Count+1; i++) {
               set s = new set();
               s.name = i;
               myset.Add(s);
           }
          
           for(int i=0;i<T.Count;i++){
               int other = T[i].Value.Key;
               int nodenum = T[i].Value.Value;
               myset[other].children.Add(myset[nodenum]);
               myset[nodenum].parent = myset[other];
           }

           return myset;
          
       }
       
    }
}