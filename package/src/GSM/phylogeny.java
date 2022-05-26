package GSM;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import fig.basic.UnorderedPair;
import goblin.Taxon;
import nuts.math.Graph;
import nuts.math.Graphs;
import nuts.util.Arbre;
import nuts.util.CollUtils;
import nuts.util.Counter;
import pty.UnrootedTree;


public class phylogeny {
	
	public static double[][] spatialtree2K(UnrootedTree urt, List<Taxon> leaves){
		int N = urt.nTaxa();
		double T = urt.totalBranchLength();

		double[][] K = new double[N][N];
		for(int i = 0; i < N; i++) {
			for(int j = i; j < N; j++) {
				if(j == i) {
					K[j][j] = 1;
				}else {
					K[i][j] = Math.exp(-3.0*urt.totalBranchLengthDistance(leaves.get(i), leaves.get(j))/T);
					K[j][i] = K[i][j];
				}
			}
		}
		return(K);
	}
	
	public static double[][] tree2K(UnrootedTree urt, List<Taxon> leaves){
		int N = urt.nTaxa();
	    int nedge =  N*2 - 3;
	    int[][] Xe = new int[N][nedge];
	    double[] mue = new double[nedge];
	    double[] sigma2e = new double[nedge];
	    double T = urt.totalBranchLength();
	    double[] e = new double[nedge];
	    
	    Counter<UnorderedPair<Set<Taxon>, Set<Taxon>>> result = new Counter<UnorderedPair<Set<Taxon>, Set<Taxon>>>();
	    List<Taxon> allLeaves = new ArrayList<Taxon>();	   
	    for(int ind = 0; ind < N; ind++) {
	    	  String ss = "leaf_"+(ind+1);
	    	  Taxon current = new Taxon(ss);
	    	  allLeaves.add(current);	    	
	    }
	    Graph<Taxon> topo = urt.getTopology();
	    Arbre<Taxon> topology = Arbre.tree2Arbre(Graphs.toTree(topo, allLeaves.get(0)));
	    Map<Arbre<Taxon>, Set<Taxon>> leavesMap = Arbre.leavesMap(topology);
	    
	    
	    int iter = 0;
	    int left = 0; int right = 0;
	    for (Arbre<Taxon> key : leavesMap.keySet())
	      if (!key.isRoot())
	      {
	        Set<Taxon> clade = leavesMap.get(key);
	        double bl = urt.branchLength(key.getContents(), key.getParent().getContents());
	        Set<Taxon> complement = CollUtils.set(allLeaves);
	        complement.removeAll(clade);
	        e[iter] = bl;
	        mue[iter] = 0.0;
    	        if(clade.size() > complement.size()) {
  	    	       left = 0; right = 1;
  	    	       for(int i = 0; i < N; i++) {
  	    	    	     if(clade.contains(allLeaves.get(i))) {
  	    	    	    	    Xe[i][iter] = 0;
  	    	    	     }else {
  	    	    	     	Xe[i][iter] = 1;
  	    	    	     }  	    	  
  	    	    	     mue[iter] += Xe[i][iter]; 	    	    	    
  	    	       }
  	    	 
  	        }else {
  	        	   left = 1; right = 0;
  	    	       for(int i = 0; i < N; i++) {
	    	    	     if(clade.contains(allLeaves.get(i))) {
	    	    	    	    Xe[i][iter] = 1;
	    	    	     }else {
	    	    	     	Xe[i][iter] = 0;
	    	    	     }	    	    	     
	    	    	     mue[iter] += Xe[i][iter];	    
	    	       }
  	    	   
  	        }
 

    		mue[iter] = mue[iter]/(N*1.0);
    		sigma2e[iter] = mue[iter]*(1-mue[iter]);
    		
	        result.setCount(new UnorderedPair<Set<Taxon>,Set<Taxon>>(complement, clade), bl);
	        iter++;
	       
	      }

		double[][] K = new double[N][N];
		for(int i = 0; i < N; i++) {
			for(int j = i; j < N; j++) {
				K[i][j] = 0.0;
				for(int iedge = 0; iedge < nedge; iedge++) {
					K[i][j] += (Xe[i][iedge]-mue[iedge])*(Xe[j][iedge]-mue[iedge])/sigma2e[iedge]*e[iedge]/T;
					
				}
				K[j][i] = K[i][j];
				
			}
			
		}
		
		
		return(K);
	}

}
