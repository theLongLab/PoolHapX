package MiscFunctions;

import java.io.*;
import java.util.*;

import PoolHap.Entrance;
import scala.concurrent.forkjoin.ThreadLocalRandom;

/**
 * @author Chen Cao 2019-08
 * Random Generate all files which could be used for PoolHapX testing.
 * Including the gold standard files, vef files, sample name file,
 * Each pool variant frequency file (In paratice, this file could be 
 * called by GATK/ SamTools). 
 * To simulate the evolution, each round n (last parameter) snps changes. 
 */  

public class RandomSimulator {
	String project_name;
	String folder_path;
	int num_pools;
	int num_total_haps;
	double ave_haps_inpool;
	int genome_len;
	int num_vars;
	int read_len;
	int outer_dist;
	int mut_each_generation;
	double ave_coverage;
	String read_type;
	int X10_fragments;
	/*	1	0	1	0	0	0	1
	 * 	1	0	0	0	0	0	1
	 * 	0	0	1	0	0	1	0
	 * 3 genomes, the 1st one is 1010001
	 */	
	int[][] sim_genome;   // num_total_haps* num_vars
	double [][] hap_freq_inpool; // num_pools* num_total_haps 
	int[] var_pos; 
	
	public int  mut(int  base)  throws IOException {
		if (base==0 ) {
			return 1;
		}else {
			return 0;
		}
	}
	
	public void GenomeSimulator() throws IOException {
		this.sim_genome = new int [this.num_total_haps][this.num_vars];
		
		for (int i = 0; i < this.sim_genome[0].length; i++) {
			Random random = new Random(); 
			int genotype  = random.nextInt(2);
			this.sim_genome[0][i] = genotype ;
		}
		
		
		for (int i = 1; i < this.sim_genome.length; i++) {
			HashSet<Integer> pos_set = new HashSet<Integer>();
			for (int j = 0; j < this.sim_genome[i].length; j++) {
				this.sim_genome[i][j] = this.sim_genome[i-1][j];
			}
			int count =0;
			while (count < this.mut_each_generation) { 
				
				Random random = new Random(); 
				int pos   = random.nextInt(this.num_vars);
				if  (!pos_set.contains(pos)){
					count++;
					sim_genome[i][pos]=mut (sim_genome[i-1][pos] );
					pos_set.add(pos);
				}
				
			}
		}	
		
		
		Boolean[] array = new Boolean[this.genome_len];
		this.var_pos = new int [this.num_vars];
		Arrays.fill(array, Boolean.FALSE);
		int count =0;
		
		while (count < this.num_vars) {
			Random random = new Random(); 
			int value   = random.nextInt(this.genome_len- 3*this.read_len-this.outer_dist)+ 
					this.read_len;
			if (array[value]==false ) {
				count++ ;
				array[value]= true ;
			}
		}
		count=0;
		for (int i = 0; i < array.length; i++) {
			if (array[i]== true) {
				this.var_pos[count]=i+1;
				count=count+1;
			}
		}
		
		String file_path= this.folder_path+"/"+ this.project_name+"/gold_standard/"+
				this.project_name+ "_haps.txt";
		FileWriter mydata = new FileWriter(file_path,false);
        PrintWriter pw = new PrintWriter(mydata);
        String pw_line ="";
        for (int i = 0; i < this.sim_genome.length; i++) {
        	pw_line="@Haplotype_"+ Integer.toString(i);
        	pw.write(pw_line+"\n");
        	pw_line="";
        	for (int j = 0; j < this.sim_genome[i].length; j++) {
        		pw_line=  pw_line+ Integer.toString(this.sim_genome[i][j]) ; 
        	}
        	pw.write(pw_line+"\n");
        }
        pw.flush();
        pw.close();
		return;	
	}
	
	public void HapFreqSimulator() throws IOException {
		double [] total_value_arr= new double [this.num_pools];
		this.hap_freq_inpool = new double [this.num_pools][this.num_total_haps];
		for (int i = 0; i < this.hap_freq_inpool.length; i++) {
			double total_value = 0.0;
			for (int j = 0; j < this.hap_freq_inpool[i].length; j++) {
				Random random = new Random(); 
				int value   = random.nextInt(10*this.num_total_haps);
				if (value< 10*(this.num_total_haps- this.ave_haps_inpool )) {
					value =0;
				}
				this.hap_freq_inpool[i][j]= value;
				total_value+= Double.valueOf(value) ;
			}
			total_value_arr[i]= total_value;
		}
		
		for (int i = 0; i < this.hap_freq_inpool.length; i++) {
			for (int j = 0; j < this.hap_freq_inpool[i].length; j++) {
				if (this.hap_freq_inpool[i][j] / total_value_arr[i] > 0.0005) {
					this.hap_freq_inpool[i][j] =  (double)Math.round(this.hap_freq_inpool[i][j] / 
							total_value_arr[i]*1000)/1000 ;
				}else {
					this.hap_freq_inpool[i][j]= 0.0;
				}
			}
		}
		
		String file_path= this.folder_path+"/"+ this.project_name+"/gold_standard/"+
				this.project_name+ "_haps.intra_freq.txt";
		FileWriter mydata = new FileWriter(file_path,false);
        PrintWriter pw = new PrintWriter(mydata); 
        String pw_line = "Hap_ID";
        for (int i = 0; i < this.num_total_haps; i++) { 
        	pw_line = pw_line+  "\t"+ "h"+Integer.toString(i);
        }
        pw.write( pw_line + "\n"); 
        for (int i = 0; i < this.hap_freq_inpool.length; i++) {
        	pw_line = this.project_name+"_p"+Integer.toString(i);
        	for (int j = 0; j < this.hap_freq_inpool[i].length; j++) {
        		pw_line =pw_line +"\t"+  Double.toString(this.hap_freq_inpool[i][j]);
        	}
        	pw.write( pw_line + "\n"); 
        }
        pw.flush();
        pw.close();
        
       
        file_path= this.folder_path+"/"+ this.project_name+"/gold_standard/"+
				this.project_name+ "_haps.inter_freq_vars.txt";
        FileWriter mydata2 = new FileWriter(file_path,false);
        PrintWriter pw2 = new PrintWriter(mydata2); 
        pw_line = "Hap_ID";
        for (int i = 0; i < this.num_total_haps; i++) { 
        	pw_line = pw_line+  "\t"+ "h"+Integer.toString(i);
        }
        pw2.write( pw_line + "\n"); 
        pw_line= "Freq";
        for (int j = 0; j < this.hap_freq_inpool[0].length; j++) {
        	double total_value = 0.0;
        	for (int i = 0; i < this.hap_freq_inpool.length; i++) {
        		total_value+=  this.hap_freq_inpool[i][j];
			}
        	pw_line=pw_line+"\t" + Double.toString(total_value/ Double.valueOf(this.num_pools)); 
		}
        pw2.write( pw_line + "\n"); 
        for (int j = 0; j < this.sim_genome[0].length; j++) {
        	pw_line = "0;"+ Integer.toString(this.var_pos[j])+";"+  Integer.toString(this.var_pos[j])
    		+";0:1";
        	for (int i = 0; i < this.sim_genome.length; i++) {
        		pw_line=pw_line + "\t"+ Integer.toString(sim_genome[i][j]) ; 
        	}
        	pw2.write( pw_line + "\n"); 
        }
        pw2.flush();
        pw2.close();
		return;
	}
	
	public void VefSimulator() throws IOException {
		if (this.read_type.equals("paired-end")) {
			new File(this.folder_path + "/"+this.project_name + "/intermediate/vef").mkdir();
			for (int i = 0; i < this.num_pools; i++) { 
				String file_path= this.folder_path+"/"+ this.project_name+"/"+
						"intermediate/vef/"+ this.project_name+ "_p"+Integer.toString(i)+".vef" ;
				
				FileWriter mydata = new FileWriter(file_path,false);
				PrintWriter pw = new PrintWriter(mydata); 
				String pw_line ="";
	//			int[][] sim_genome;   // num_total_haps* num_vars
	//			double [][] hap_freq_inpool; // num_pools* num_total_haps
				int [] hap_count = new int [this.num_total_haps ];
				int roulette_total =0;
				HashMap<Integer, Integer> hap_dict = new HashMap<Integer, Integer>(); 
				for (int j = 0; j < this.num_total_haps; j++) { 
					int count = (int) (hap_freq_inpool[i][j]*1000);
					hap_count[j]= count ;	
					for (int k = 0; k < count ; k++) {
						hap_dict.put(roulette_total, j);
						roulette_total++;
					}
				}
				
				int total_read = (int) ((this.ave_coverage* this.genome_len)/ (this.read_len*2 ));
				for (int t = 0; t < total_read; t++) { 
					Random random = new Random(); 
	//				System.out.println(roulette_total);
					int value   = random.nextInt( roulette_total ) ;
					int hap = hap_dict.get(value);
					Random random_pos = new Random(); 
					int pos =  random_pos.nextInt(this.genome_len ) ;
	//				@Haplotype_0_4618_4869_0_1_0_0_0:0:0_0:0:0_1c/1
					int start_1=pos;
					int end_1=pos+this.read_len;
					int start_2=pos+ this.read_len+ this.outer_dist;
					int end_2=pos+ 2*this.read_len+ this.outer_dist;
					pw_line="@Haplotype_"+ Integer.toString(hap)+"_"+ Integer.toString(t)+"\t";
					boolean flag= false;
					HashMap<Integer, Integer> var_dict = new HashMap<Integer, Integer>();
					for (int j = 0; j < this.num_vars; j++) { 
						var_dict.put(var_pos[j],sim_genome[hap][j] );
					}
					for (int j = start_1; j < end_1 ; j++) {
						if (var_dict.containsKey(j)) {
							pw_line=pw_line+ Integer.toString(j)+"="+var_dict.get(j)+";";
							flag=true;
						}
					}
					for (int j = start_2; j < end_2 ; j++) {
						if (var_dict.containsKey(j)) {
							pw_line=pw_line+ Integer.toString(j)+"="+var_dict.get(j)+";";
							flag=true;
						}
					}
					if (flag==true ) {
						pw_line=pw_line+ "\t"+ "//\t"+ Integer.toString(start_1)+"\t"+ 
								Integer.toString(end_1) +"\t"+
								Integer.toString(start_2)+"\t"+Integer.toString(end_2);
						pw.write( pw_line + "\n"); 
					}
				}
				pw.flush();
		        pw.close();
				
			}
		}else if  (this.read_type.equals("10X")) {
			this.X10_fragments= 10;
			new File(this.folder_path + "/"+this.project_name + "/intermediate/vef").mkdir();
			for (int i = 0; i < this.num_pools; i++) { 
				String file_path= this.folder_path+"/"+ this.project_name+"/"+
						"intermediate/vef/"+ this.project_name+ "_p"+Integer.toString(i)+".vef" ;
				
				FileWriter mydata = new FileWriter(file_path,false);
				PrintWriter pw = new PrintWriter(mydata); 
				String pw_line ="";
				//			int[][] sim_genome;   // num_total_haps* num_vars
				//			double [][] hap_freq_inpool; // num_pools* num_total_haps
				int [] hap_count = new int [this.num_total_haps ];
				int roulette_total =0;
				HashMap<Integer, Integer> hap_dict = new HashMap<Integer, Integer>(); 
				for (int j = 0; j < this.num_total_haps; j++) { 
					int count = (int) (hap_freq_inpool[i][j]*1000);
					hap_count[j]= count ;	
					for (int k = 0; k < count ; k++) {
						hap_dict.put(roulette_total, j);
						roulette_total++;
					}
				}
				int total_read = (int) ((this.ave_coverage* this.genome_len)/ 
						(this.read_len*this.X10_fragments ));
				for (int t = 0; t < total_read; t++) { 
					Random random = new Random(); 
					int value   = random.nextInt( roulette_total ) ;
					int hap = hap_dict.get(value);
					Random random_pos = new Random(); 
					int pos =  random_pos.nextInt(this.genome_len ) ;
					pw_line="@Haplotype_"+ Integer.toString(hap)+"_"+ Integer.toString(t)+"\t";
					int [] start_vec = new int [ this.X10_fragments];
					int [] end_vec = new int [ this.X10_fragments];
					start_vec[0]= pos;
					end_vec[0]= pos+this.read_len;
					for (int p=1; p< start_vec.length;p++) {
						start_vec[p]= end_vec[p-1]+ this.outer_dist;
						end_vec[p]= start_vec[p]+ this.read_len;
					}
					boolean flag= false;
					HashMap<Integer, Integer> var_dict = new HashMap<Integer, Integer>();
					for (int j = 0; j < this.num_vars; j++) { 
						var_dict.put(var_pos[j],sim_genome[hap][j] );
					}
					for (int f = 0; f < this.X10_fragments ; f++) {
						for (int j=start_vec[f]; j<end_vec[f];j++) {
							if (var_dict.containsKey(j)) {
								pw_line=pw_line+ Integer.toString(j)+"="+var_dict.get(j)+";";
								flag=true;
							}
						}
					}
					
					if (flag==true ) {
						pw_line=pw_line+ "\t"+ "//\t"+ Integer.toString(start_vec[0])+"\t"+ 
								Integer.toString(end_vec[0]) +"\t"+
								Integer.toString(start_vec[this.X10_fragments-1])+"\t"+
									Integer.toString(end_vec[this.X10_fragments-1]);
						pw.write( pw_line + "\n"); 
					}
				}
				
				
				pw.flush();
		        pw.close();
			}
			
		}
	}
	
	public void VarFreqSimulator() throws IOException {
		String file_path= this.folder_path+"/"+ this.project_name+"/"+"intermediate/"
				+ this.project_name+ "_sample_names.txt";
		FileWriter mydata0 = new FileWriter(file_path,false);
        PrintWriter pw0 = new PrintWriter(mydata0); 
        for (int i = 0; i < this.num_pools; i++) { 
        	pw0.write( this.project_name +"_p"+ Integer.toString(i ) + "\n");  
        }
        pw0.flush();
        pw0.close();
        
        file_path= this.folder_path+"/"+ this.project_name+"/"+"intermediate/" +
				this.project_name+ "_vars.intra_freq.txt";
		FileWriter mydata = new FileWriter(file_path,false);
		PrintWriter pw = new PrintWriter(mydata); 
		String pw_line="Pool_ID";
		for (int i = 0; i < this.num_pools; i++) { 
			pw_line=pw_line+ "\t"+ this.project_name +"p"+ Integer.toString(i ) ;
		}
		pw.write(pw_line+"\n");
		for (int i = 0; i < this.num_vars; i++) {
			pw_line = "0;"+ Integer.toString(this.var_pos[i])+";"+  Integer.toString(this.var_pos[i])
    		+";0:1";
			for (int j = 0; j < this.num_pools; j++) {
				double freq = 0.0;
//				int[][] sim_genome;   // num_total_haps* num_vars
//				double [][] hap_freq_inpool; // num_pools* num_total_haps
				for (int k = 0; k < this.num_total_haps; k++) {
					freq+= hap_freq_inpool[j][k]* sim_genome [k][i];
				}
				freq= (double)Math.round(freq *1000)/1000;
				pw_line=pw_line + "\t"+ Double.toString(freq )  ;
			}
			pw.write(pw_line+"\n");
		}
        
		pw.flush();
        pw.close();
	}
	
	public void LassoSimulator() throws IOException {
//		ArrayList<String >  random_haps = new ArrayList<String>();
//		HashSet<Integer> haps_set = new HashSet<Integer>();
		int[][] random_genome = new int [this.num_total_haps* (this.num_vars-1)]
				[this.num_vars];
		int index =0;
		for (int i = 0; i < this.num_total_haps; i++) { 
			int [] tmp_genome2 = new int [this.num_vars];
			for (int k=0;k<this.sim_genome[i].length;k++ ) {
				tmp_genome2[k]= this.sim_genome[i][k];
			}
			for (int k=0;k<tmp_genome2.length;k++ ) {
				random_genome[index][k] = tmp_genome2[k];
			}
			index ++;
			for (int j = 0; j < (this.num_vars-2); j++) {
				int [] tmp_genome = new int [this.num_vars];
				for (int k=0;k<this.sim_genome[i].length;k++ ) {
					tmp_genome[k]= this.sim_genome[i][k];
				}
				tmp_genome[j] =this.mut(this.sim_genome[i][j]); 
				tmp_genome[j+1] =this.mut(this.sim_genome[i][j+1]); 
				tmp_genome[j+2] =this.mut(this.sim_genome[i][j+1]); 
				
//				for (int k=0;k<this.sim_genome[i].length;k++ ) {
//					Random random = new Random(); 
//					int genotype  = random.nextInt(2);
//					tmp_genome[k]= genotype;
//				}
				for (int k=0;k<tmp_genome.length;k++ ) {
					random_genome[index][k] = tmp_genome[k];
				}
				index ++;
			}
		}
		
		
		String file_path2= this.folder_path+"/"+ this.project_name+"/"+
				"output/"+ this.project_name+ "_randome_genomes.txt" ;
		FileWriter mydata2 = new FileWriter(file_path2,false);
		PrintWriter pw2 = new PrintWriter(mydata2); 
		String pw2_line ="";
		for (int i = 0; i < random_genome.length; i++) { 
			pw2_line="";
			for (int j = 0; j < random_genome[i].length; j++) { 
				pw2_line=pw2_line +random_genome[i][j];  
			}
			pw2.write(pw2_line+"\n");
		}
		
		
		pw2.flush();
        pw2.close();
        
		
		new File(this.folder_path + "/"+this.project_name + "/intermediate/lasso/").mkdir();
		for (int i = 0; i < this.num_pools; i++) { 
//		for (int i = 4; i < 5; i++) { 
			String file_path= this.folder_path+"/"+ this.project_name+"/"+
					"intermediate/lasso/"+ this.project_name+ "_p"+Integer.toString(i)+".lasso_in" ;
			FileWriter mydata = new FileWriter(file_path,false);
			PrintWriter pw = new PrintWriter(mydata); 
			String pw_line ="1.0";
			for (int j = 0; j < random_genome.length; j++) { 
				pw_line=pw_line +" "+ Integer.toString(j+1)+":1.0";
			}
			pw.write( pw_line + "\n");  
			pw_line="";
			
			
			
//			for (int j = 0; j < this.num_vars; j++) { 
//				System.out.print(this.sim_genome[0][j]); 
//			}
//			System.out.println();
			for (int j = 0; j < this.num_vars; j++) { 
				double freq =0;
				for (int k=0; k< this.hap_freq_inpool[i].length; k++ ) {
					if (this.hap_freq_inpool[i][k]> 0.001) {
						freq+= this.hap_freq_inpool[i][k]* this.sim_genome[k][j];
					}
				}
				pw_line= Double.toString(freq);
				for (int k=0; k< random_genome.length; k++ ) {
					if (random_genome[k][j]==1) {
						pw_line=pw_line+" "+ Integer.toString(k+1)+":1.0";
					}else {
						pw_line=pw_line+" "+ Integer.toString(k+1)+":0.0";
					}
				}
				pw.write( pw_line + "\n"); 
			}
			
			
			for (int j = 0; j < (this.num_vars-1 ); j++) { 
				for (int k = j+1; k < this.num_vars; k++) { 
					double freq =0; 
					for (int h=0; h< this.hap_freq_inpool[i].length; h++ ) {
						if (this.hap_freq_inpool[i][h]> 0.001) {
							freq+= this.hap_freq_inpool[i][h]* this.sim_genome[h][j]
									*this.sim_genome[h][k] ;
						}
					}
					if (freq >-0.0001) {
						pw_line= Integer.toString(j)+"\t"+ Integer.toString(k)+"\t"+
								Double.toString(freq);
						pw_line= Double.toString(freq);
						for (int g=0; g< random_genome.length; g++ ) {
							if( (random_genome[g][j]==1) &&  (random_genome[g][k]==1) ){
								pw_line=pw_line+" "+ Integer.toString(g+1)+":1.0";
							}else {
								pw_line=pw_line+" "+ Integer.toString(g+1)+":0.0";
							}
						}
						pw.write( pw_line + "\n"); 
					}
				}
			}
			
			pw.flush();
	        pw.close();
			
		}
		String file_path= this.folder_path+"/"+ this.project_name+"/"+
				"output/"+ this.project_name+ "_gc.inter_freq_vars.txt" ;
		FileWriter mydata = new FileWriter(file_path,false);
		PrintWriter pw = new PrintWriter(mydata); 
		String pw_line ="Hap_ID";
		for (int j = 0; j < random_genome.length; j++) { 
			pw_line=pw_line+"\t"+"h"+ Integer.toHexString(j);  ;
		}
		pw.write(pw_line+"\n"); 
		pw_line ="Freq";
		for (int j = 0; j < random_genome.length; j++) { 
			pw_line=pw_line+"\t"+ Double.toString(1/(double)random_genome.length );
		}
		pw.write(pw_line+"\n"); 
		for (int j = 0; j < this.num_vars; j++) { 
			pw_line= "0;"+ Integer.toString(this.var_pos[j])+";"+ Integer.toString(this.var_pos[j])+
					";0:1";
			for (int i=0;i<random_genome.length; i++ ) {
				pw_line =pw_line+"\t"+ Integer.toString(random_genome[i][j]); 
			}
			pw.write(pw_line+"\n"); 
			
		}
//		0;737;737;0:1
		
//		pw.write(pw_line+"\n"); 
		pw.flush();
        pw.close();
	}
	
	
	public void FastaSimulator() throws IOException {

		String fasta_file= this.folder_path+"/"+ this.project_name+"/fasta/"+"bac.fa";
		String genome="";
		BufferedReader bufferedreader = new BufferedReader(new FileReader(fasta_file));
        String line = "";
        while ((line = bufferedreader.readLine()) != null) {
        	line=line.replace("\n", "").replace("\r", "");
        	if (!line.startsWith(">")) {
        		genome=genome+ line;
        	}
        }
        bufferedreader.close();
        for (int i=0;i< this.num_total_haps;i++) {
        	BufferedWriter bw = new BufferedWriter(new FileWriter(this.folder_path+"/"+ 
        			this.project_name+"/fasta/"+"bac."+ Integer.toString(i)+".fa"));
        	bw.write(">>NC_014932.1_h"+Integer.toString(i)+ " Bartonella clarridgeiae strain 73, complete genome\n");
        	String tmp_genome =genome;
        	for (int j=0;j< sim_genome[i].length;j++) {
        		if (sim_genome[i][j]==1) {
        			int pos= var_pos[j];
        			pos=pos-1;
        			if (tmp_genome.substring(pos,pos+1).equals("A") ) {
        				tmp_genome= tmp_genome.substring(0, pos)+"T"+ tmp_genome.substring(pos+1);
        			}else {
        				tmp_genome= tmp_genome.substring(0, pos)+"A"+ tmp_genome.substring(pos+1);
        			}
        		}
        	}
        	bw.write( tmp_genome);
        	bw.close();
        	
        }
        
        for (int i=0;i< this.num_pools;i++) {
        	System.out.println("mkdir /home/chencao/Desktop/bac/fastq/bac_p"+Integer.toString(i));
        	System.out.println("cd /home/chencao/Desktop/bac/fastq/bac_p"+Integer.toString(i));
//        	hap_freq_inpool
        	for (int j=0; j<hap_freq_inpool[i].length; j++ ) {
        		if (hap_freq_inpool[i][j]>0.001 ) {
        			String tmp = "dwgsim /home/chencao/Desktop/bac/fasta/bac."+Integer.toString(j) +".fa"+
        					" -e 0.01 -E 0.01 -C "+ Double.toString( hap_freq_inpool[i][j]*100) + " -1 150 -2 150 -r 0 "
        							+ " bac.h"+ Integer.toString(j);
        			System.out.println(tmp);
        		}
        	}
        }
	}
	
	
	public RandomSimulator(String[] args) throws IOException {
		
 

	
		this.folder_path= args[0];
		this.project_name= args[1];
		this.num_pools = Integer.parseInt(args[2]);
		this.num_total_haps = Integer.parseInt(args[3]);
		this.ave_haps_inpool= Double.parseDouble( args[4]);
		this.genome_len= Integer.parseInt(args[5]);
		this.num_vars = Integer.parseInt(args[6]);
		this.read_len = Integer.parseInt(args[7]);
		this.outer_dist = Integer.parseInt(args[8]);
		this.ave_coverage= Double.parseDouble( args[9]);
		this.mut_each_generation = Integer.parseInt(args[10]);
		this.read_type= args[11];
		new File(this.folder_path + "/"+this.project_name ).mkdir();
		new File(this.folder_path + "/"+this.project_name + "/"+ "/gold_standard/").mkdir();
		new File(this.folder_path + "/"+this.project_name + "/"+ "/input/").mkdir();
		new File(this.folder_path + "/"+this.project_name + "/"+ "/output/").mkdir();
		new File(this.folder_path + "/"+this.project_name + "/"+ "/intermediate/").mkdir();
		GenomeSimulator (); 
		System.out.println(" Genome Simulation Finished!");
		HapFreqSimulator();
		System.out.println(" Haplotype Frequency for Each Pool Simulation Finished!");
		VarFreqSimulator();
		System.out.println(" Variants Frequency for Each Haplotype Simulation Finished!");
		VefSimulator();
		System.out.println(" Vef File Simulation Finished!");
		
//		FastaSimulator();
//		System.out.println(" Fasta Simulation Finished!");
//		LassoSimulator();
//		System.out.println(" Lasso Input File Simulation Finished!");
		

		
	}
	
	

	
	public static  void main(String[] args) throws IOException, InterruptedException {
		int [] hap_count = new int [51]; 
		for (int i=0;i<1100;i++) {
			int x= ThreadLocalRandom.current().nextInt(0, 50);
			hap_count[x]++;
//			System.out.println(ThreadLocalRandom.current().nextInt(0, 50));
		}
		for (int i=0;i<51;i++) {
			System.out.println(i+"\t"+hap_count[i]);
		}
		
//		String x="012345";
//		String y="67890123";
//		String z="4567890";
//		String tmp =x+y;
//		System.out.println(tmp.substring(0, 4) + z+ tmp.substring(4+7)) ;
//		for (int h=0;h< 100;h++) {
//			Random seedRandom = new Random( h);
//	        int currPool = seedRandom.nextInt( 20);
//	        System.out.println(currPool);
//		}
		
//		Process ps = Runtime.getRuntime().exec("pwd ", null);
//		Process ps = Runtime.getRuntime().exec("du -sh ", null);
//		
//        ps.waitFor();  
//        BufferedReader br = new BufferedReader(new InputStreamReader(ps.getInputStream()));  
//        StringBuffer sb = new StringBuffer();  
//        String line;  
//        while ((line = br.readLine()) != null) {  
//            sb.append(line).append("\n");  
//        }  
//    	String result = sb.toString();  
//    	System.out.println(result);  
//    	br.close();
		
		
		
		//		System.exit(0);
		// /home/chencao/Desktop  sim001	10	20	15	5000	50	150	50	100		5
		//parameter 0: folder path
		//parameter 1: project name
		//parameter 2: number of pools (or generations)
		//parameter 3: number of total haplotypes
		//parameter 4: average number of haplotypes in each pool
		//parameter 5: genome length 
		//parameter 6: number of variants
		//parameter	7: read length
		//parameter 8: outer distance between the two ends for pairs
		//parameter 9: average coverage
		//parameter 10: number of mutations each generation
		//parameter 11: 10X or paired-end
		RandomSimulator rs = new RandomSimulator(args);
	}
	

}
