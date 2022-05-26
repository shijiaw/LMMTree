package GSM;

import java.io.File;
import java.io.PrintWriter;
import java.io.*;  
import java.util.Scanner;  
import java.util.ArrayList;
import java.util.List;

import fig.basic.IOUtils;
import goblin.Taxon;
import nuts.io.CSV;
import nuts.io.IO;
import pty.RootedTree;
import pty.UnrootedTree;

public class gsm {
	public static void main(String[] args) throws Exception   {
	
	    Scanner sc = new Scanner(new File("../M.txt"));  
        sc.useDelimiter(",");   //sets the delimiter pattern  
        double M = Double.parseDouble(sc.next());  
        sc.close(); 
		
		long time = System.currentTimeMillis();

		for(int index = 0; index < M; index++) {
    		String inputTree = "../mb_tree/tree_"+ (index+1) + ".newick";
            String OutputFileName = "../mb_tree/K_"+ (index+1) + ".csv";
    		
			PrintWriter KOut = IOUtils.openOutEasy(new File(OutputFileName));

			File dataFile = new File(inputTree);
			String s = IO.f2s(dataFile);
			System.out.println(s);
			UnrootedTree goldut = UnrootedTree.fromRooted(RootedTree.Util.fromNewickString(IO.f2s(dataFile)));
			List<Taxon> leaves = goldut.leaves();
			
			String st = dataFile.getName();
			double[][] K = phylogeny.tree2K(goldut, leaves);
			
			for(int i = 0; i < leaves.size(); i++) {
				List<Double> Kind = new ArrayList<>();
				for(int ind = 0; ind < leaves.size(); ind++) {
					Kind.add(K[i][ind]);
				}
				KOut.println(CSV.body(Kind));
			}		
			
			KOut.close();
		}
		time = System.currentTimeMillis() - time;
		System.out.println(time);
	}
}
