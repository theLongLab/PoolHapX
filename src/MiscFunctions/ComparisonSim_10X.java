package MiscFunctions;

import java.io.*;
import java.util.*;

import PoolHap.Entrance;

public class ComparisonSim_10X {
	
//	int[] pools = new  int [] {50, 200};
	int[] pools = new  int [] {25, 50};
//	double[] mut_rates = new  double [] { 1e-6, 1e-7, 1e-8, 1e-9};
//	int[] depths = new  int  [] { 250, 500, 1000, 2000};
	int[] depths = new  int  [] {100, 250 };
	int [] freq_cutoff =new  int  [] {0, 2, 4 };
//	int[] num_haps = new  int  [] { 15, 30};
	
	String prefix_folder= "";
	String slim_script;  
// 
	int genome_len =  1067971 ;
	String genome_path ="/export/qlong/chencao/Work/poolhapx/slim/ref/plasmodium/plasmodium.fa"; 
	String species= "10X"; //virus, bacteria, metagenomics, human
	
	
	public ComparisonSim_10X(String parm ) throws IOException {
		
		this.prefix_folder = parm+"/";
		new File(this.prefix_folder + "/cmd/").mkdir();
		
		BufferedWriter bw_sh = new BufferedWriter(new FileWriter(this.prefix_folder + "/cmd/cmd.sh"));
		
		for (int f =0; f< this.freq_cutoff.length; f++) {
			for (int i =0; i< this.pools.length; i++) {
				for (int j =0; j< this.depths.length; j++) {
					bw_sh.write("sbatch freq_" +Integer.toString(this.freq_cutoff[f])+ 
							"_pool_"+Integer.toString(this.pools[i])+ "_dep_"+ 
							Integer.toString(this.depths[j]) +".1.cmd\n"  );
				}
			}
		}
		
		bw_sh.close();
		
		BufferedWriter bw_sh2 = new BufferedWriter(new FileWriter(this.prefix_folder + "/cmd/cmd2.sh"));
		
		for (int f =0; f< this.freq_cutoff.length; f++) {
			for (int i =0; i< this.pools.length; i++) {
				for (int j =0; j< this.depths.length; j++) {
					bw_sh2.write("sbatch freq_" +Integer.toString(this.freq_cutoff[f])+ 
							"_pool_"+Integer.toString(this.pools[i])+ "_dep_"+ 
							Integer.toString(this.depths[j]) +".2.cmd\n"  );
				}
			}
		}
		
		bw_sh2.close();
		
		
		BufferedWriter bw_sh3 = new BufferedWriter(new FileWriter(this.prefix_folder + "/cmd/cmd3.sh"));
		
		for (int f =0; f< this.freq_cutoff.length; f++) {
			for (int i =0; i< this.pools.length; i++) {
				for (int j =0; j< this.depths.length; j++) {
					bw_sh3.write("sbatch freq_" +Integer.toString(this.freq_cutoff[f])+ 
							"_pool_"+Integer.toString(this.pools[i])+ "_dep_"+ 
							Integer.toString(this.depths[j]) +".3.cmd\n"  );
				}
			}
		}
		
		bw_sh3.close();
		
		
		
		for (int f =0; f< this.freq_cutoff.length; f++) {
			for (int i =0; i< this.pools.length; i++) {
				for (int j =0; j< this.depths.length; j++) {
					new File(this.prefix_folder + "/freq_"+ Integer.toString(this.freq_cutoff[f])+ 
							"_pool_"+ Integer.toString(this.pools[i])
							+ "_dep_"+ Integer.toString(this.depths[j])).mkdir();
					
					BufferedWriter bw1 = new BufferedWriter(new FileWriter(this.prefix_folder + "/cmd/freq_"
							+ Integer.toString(this.freq_cutoff[f])
							+ "_pool_"+ 
							Integer.toString(this.pools[i])+ "_dep_"+ 
							Integer.toString(this.depths[j]) +".PoolSimulator.properties"));
							
					bw1.write("Input_Dir = "+ this.prefix_folder + "/freq_"+ Integer.toString(this.freq_cutoff[f])
							+ "_pool_"+ Integer.toString(this.pools[i])
					+ "_dep_"+ Integer.toString(this.depths[j])+ "/input\n");
					bw1.write("Intermediate_Dir = "+ this.prefix_folder + "/freq_"+ Integer.toString(this.freq_cutoff[f])
					+ "_pool_"+ Integer.toString(this.pools[i])
					+ "_dep_"+ Integer.toString(this.depths[j])+ "/intermediate\n");
					bw1.write("Gold-Standard_Dir = "+ this.prefix_folder + "/freq_"+ Integer.toString(this.freq_cutoff[f])
					+ "_pool_"+ Integer.toString(this.pools[i])
					+ "_dep_"+ Integer.toString(this.depths[j])+ "/gold_standard\n");
					bw1.write("Slim_Output_Path = "+ this.prefix_folder + "/freq_"+ Integer.toString(this.freq_cutoff[f])
					+ "_pool_"+ Integer.toString(this.pools[i])
					+ "_dep_"+ Integer.toString(this.depths[j])+ "/gold_standard\n");

	//				Slim_Model = panmictic_haploid
					bw1.write("Slim_Model = island_haploid\n");
					bw1.write("Proj_Name = freq_"+  Integer.toString(this.freq_cutoff[f])+
							"_pool_" + Integer.toString(this.pools[i]) 
							+ "_dep_"+ Integer.toString(this.depths[j])+"\n");
					bw1.write("Is_Single_Population = false\n");
					bw1.write("Is_Ms_Output = false\n");
					bw1.write("DWGSIM = /export/home/jhe/download/DWGSIM-master/dwgsim\n");
					bw1.write("Num_Pools = " + Integer.toString(this.pools[i])+ "\n");
					bw1.write("Ref_Seq_Len = "+ Integer.toString(this.genome_len)+"\n");
					bw1.write("Reference_Seq= " + this.genome_path+"\n");
					bw1.write("Is_Perfect = false\n");
					bw1.write("Error_Rate_Per_Base = 0.001 \n");
					bw1.write("Hap_Freq_Cutoff = "+ 
					Double.toString(Double.valueOf(this.freq_cutoff[f])*0.005)+"   \n");
					String cov= Integer.toString(this.depths[j]);
					if (cov.equals("250")) {
						cov="300";
					}
					bw1.write("Coverage = "+ cov+  "\n");
					bw1.write("Read_Len = 150\n");
					bw1.write("Outer_Dist = 400\n");
					bw1.write("Weak_Length = 400\n");
					bw1.write("Num_Haps_Pool= 30\n");
//					bw1.write("Num_Pools= "+Integer.toString(this.pools[i])+   "\n");
					bw1.close();
					
					BufferedWriter bw = new BufferedWriter(new FileWriter(this.prefix_folder + "/cmd/freq_"+
							Integer.toString(this.freq_cutoff[f])+"_pool_"+  
							Integer.toString(this.pools[i])
							+ "_dep_"+ Integer.toString(this.depths[j]) +".1.cmd"));
					
					bw.write("#!/bin/bash\n");
					bw.write("#SBATCH --job-name=f_"+Integer.toString(this.freq_cutoff[f])+
							"_p_"+ Integer.toString(this.pools[i])
							+ "_d_"+ Integer.toString(this.depths[j]) +"\n");
					bw.write("#SBATCH --workdir="+ this.prefix_folder + "/freq_"+
							Integer.toString(this.freq_cutoff[f])+ "_pool_"+ Integer.toString(this.pools[i])
					+ "_dep_"+ Integer.toString(this.depths[j]) +"\n");
					bw.write("#SBATCH --error="+this.prefix_folder+ "/cmd/f_"+ Integer.toString(this.freq_cutoff[f])+
							"_p_"+Integer.toString(this.pools[i])
					+ "_d_"+ Integer.toString(this.depths[j])+".1.error\n" );
					bw.write("#SBATCH --output="+this.prefix_folder+ "/cmd/f_"+ Integer.toString(this.freq_cutoff[f])
							+ "_p_"+Integer.toString(this.pools[i])
					+ "_d_"+ Integer.toString(this.depths[j])+".1.output\n" );
					
					bw.write("#SBATCH --ntasks=1\n");
					bw.write("#SBATCH --cpus-per-task=6\n");
					bw.write("#SBATCH --time=99-00:00:00\n");
					bw.write("#SBATCH --nodes=1\n");
					
					
	//				java=/home/jingni.he1/download/java_jdk_8u201/jdk1.8.0_201/bin/java
	//				slim=/home/jingni.he1/download/SliM/build/slim
					bw.write("java=/export/home/jhe/download/java_jdk_8u201/jdk1.8.0_201/bin/java\n");
					bw.write("slim=/export/home/jhe/download/SliM/build/slim\n");
					bw.write("poolhapx=/export/qlong/PoolHapX/PoolHapX.jar\n");
					if (!species.equals("bacteria")) {
//						bw.write("poolsim=/export/home/jhe/project/Viral_reconstruction/SLiM/"
//							+ "programs/PoolSimulator_SLiM.jar\n");
						bw.write("poolsim=/export/qlong/PoolHapX/PoolSimulator_SLiM_10X.jar\n");
					} else{
						bw.write("poolsim=/export/qlong/PoolHapX/PoolSimulator_SLiM_bacteria.jar\n");
					}
					
					bw.write("bwa=/export/home/jhe/download/bwa-0.7.17/bwa\n");
					bw.write("ref="+ this.genome_path+"\n");
					bw.write("gatk=/export/qlong/chencao/Work/poolhapx/software/gatk-4.0.0.0/"
							+ "gatk-package-4.0.0.0-local.jar\n");
					
	//				/home/jingni.he1/project/Viral_reconstruction/SLiM/programs/PoolSimulator_SLiM.jar	
	//$java -jar /home/jingni.he1/project/Viral_reconstruction/SLiM/programs/
	//				PoolSimulator_SLiM.jar "$1"/input/PoolSimulator.properties
					
					
					bw.write("rm -rf  "+ this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])+
							"_pool_"+ Integer.toString(this.pools[i])
					+ "_dep_"+ Integer.toString(this.depths[j])+ "/input/fastq\n");
					
					bw.write("rm -rf  "+ this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])+
							"_pool_"+ Integer.toString(this.pools[i])
					+ "_dep_"+ Integer.toString(this.depths[j])+ "/input/fasta\n");
					
					bw.write("rm -rf  "+ this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])+
							"_pool_"+ Integer.toString(this.pools[i])
					+ "_dep_"+ Integer.toString(this.depths[j])+ "/input/bam\n");
					
					
					bw.write("rm -rf  "+ this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])+
							"_pool_"+ Integer.toString(this.pools[i])
					+ "_dep_"+ Integer.toString(this.depths[j])+ "/input/sam\n");
					
					bw.write("rm -rf  "+ this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])+
							"_pool_"+ Integer.toString(this.pools[i])
					+ "_dep_"+ Integer.toString(this.depths[j])+ "/input/vcf\n");
					
					bw.write("rm -rf  "+ this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])+
							"_pool_"+ Integer.toString(this.pools[i])
					+ "_dep_"+ Integer.toString(this.depths[j])+ "/gold_standard\n");
					
					
					bw.write("rm -rf  "+ this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])+
							"_pool_"+ Integer.toString(this.pools[i])
					+ "_dep_"+ Integer.toString(this.depths[j])+ "/output\n");
					
					bw.write("rm -rf  "+ this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])+
							"_pool_"+ Integer.toString(this.pools[i])
					+ "_dep_"+ Integer.toString(this.depths[j])+ "/intermediate\n");
					
				
					
					
					bw.write("mkdir "+ this.prefix_folder + "/freq_"+ Integer.toString(this.freq_cutoff[f])+ 
					"_pool_"+ Integer.toString(this.pools[i])
					+ "_dep_"+ Integer.toString(this.depths[j])+ "/input\n");
					
					bw.write("mkdir "+ this.prefix_folder + "/freq_"+ Integer.toString(this.freq_cutoff[f])+ 
							"_pool_"+ Integer.toString(this.pools[i])
							+ "_dep_"+ Integer.toString(this.depths[j])+ "/gold_standard\n");
					
					bw.write("mkdir "+ this.prefix_folder + "/freq_"+ Integer.toString(this.freq_cutoff[f])+ 
							"_pool_"+ Integer.toString(this.pools[i])
							+ "_dep_"+ Integer.toString(this.depths[j])+ "/output\n");
					
					bw.write("mkdir "+ this.prefix_folder + "/freq_"+ Integer.toString(this.freq_cutoff[f])+ 
							"_pool_"+ Integer.toString(this.pools[i])
							+ "_dep_"+ Integer.toString(this.depths[j])+ "/intermediate\n");
					
					
	
					
					bw.write("cp  "+ this.prefix_folder + "island_haploid.out "
							+ this.prefix_folder + "/freq_"+ Integer.toString(this.freq_cutoff[f])+ "_pool_"+ Integer.toString(this.pools[i])
					+ "_dep_"+ Integer.toString(this.depths[j])+ "/gold_standard/freq_"+ Integer.toString(this.freq_cutoff[f])+
					"_pool_"+ Integer.toString(this.pools[i])+ "_dep_"+ Integer.toString(this.depths[j])
							+ "_island_haploid.out\n");
					
	// Step 1: Generate coalescence-simulated haplotypes, distribute to each of the pools, and simulate reads for each pool.
					
					bw.write("$java -jar $poolsim "+ 	this.prefix_folder + "/cmd/freq_"+ Integer.toString(this.freq_cutoff[f])
							+ "_pool_"+ 
					Integer.toString(this.pools[i])+ "_dep_"+ 
							Integer.toString(this.depths[j]) +".PoolSimulator.properties\n");
					bw.close();
					
					
					
//---------------------------------------------------------------------------------------------------------------------
					BufferedWriter bw2 = new BufferedWriter(new FileWriter(this.prefix_folder + "/cmd/freq_"+
							Integer.toString(this.freq_cutoff[f])+"_pool_"+  
							Integer.toString(this.pools[i])
							+ "_dep_"+ Integer.toString(this.depths[j]) +".2.cmd"));
					
					bw2.write("#!/bin/bash\n");
					bw2.write("#SBATCH --job-name=f_"+Integer.toString(this.freq_cutoff[f])+
							"_p_"+ Integer.toString(this.pools[i])
							+ "_d_"+ Integer.toString(this.depths[j]) +"\n");
					bw2.write("#SBATCH --workdir="+ this.prefix_folder + "/freq_"+
							Integer.toString(this.freq_cutoff[f])+ "_pool_"+ Integer.toString(this.pools[i])
					+ "_dep_"+ Integer.toString(this.depths[j]) +"\n");
					bw2.write("#SBATCH --error="+this.prefix_folder+ "/cmd/f_"+ Integer.toString(this.freq_cutoff[f])+
							"_p_"+Integer.toString(this.pools[i])
					+ "_d_"+ Integer.toString(this.depths[j])+".2.error\n" );
					bw2.write("#SBATCH --output="+this.prefix_folder+ "/cmd/f_"+ Integer.toString(this.freq_cutoff[f])
							+ "_p_"+Integer.toString(this.pools[i])
					+ "_d_"+ Integer.toString(this.depths[j])+".2.output\n" );
					
					bw2.write("#SBATCH --ntasks=1\n");
					bw2.write("#SBATCH --cpus-per-task=24\n");
					bw2.write("#SBATCH --time=99-00:00:00\n");
					bw2.write("#SBATCH --nodes=1\n");
					
					bw2.write("java=/export/home/jhe/download/java_jdk_8u201/jdk1.8.0_201/bin/java\n");
				
					bw2.write("ref="+ this.genome_path+"\n");
					bw2.write("gatk=/export/qlong/chencao/Work/poolhapx/software/gatk-4.0.0.0/"
							+ "gatk-package-4.0.0.0-local.jar\n");
					
					bw2.write("java=/export/home/jhe/download/java_jdk_8u201/jdk1.8.0_201/bin/java\n");
					
					bw2.write("java=/export/home/jhe/download/java_jdk_8u201/jdk1.8.0_201/bin/java\n");
					
					for (int p=0;p < this.pools[i];p ++) {
						bw2.write("prefix="+ this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])+ 
								"_pool_"+ Integer.toString(this.pools[i])+ 
								"_dep_"+ Integer.toString(this.depths[j])+"/input/fastq/freq_"+ Integer.toString(this.freq_cutoff[f])
								+ "_pool_"+ Integer.toString(this.pools[i])+ 
								"_dep_"+ Integer.toString(this.depths[j])+
								"_p"+ Integer.toString(p) +"\n");
						bw2.write("mkdir  $prefix\n");
						
						bw2.write("prefix_fa="+ this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])+ 
								"_pool_"+ Integer.toString(this.pools[i])+ 
								"_dep_"+ Integer.toString(this.depths[j])+"/input/fasta/freq_"+ Integer.toString(this.freq_cutoff[f])
								+ "_pool_"+ Integer.toString(this.pools[i])+ 
								"_dep_"+ Integer.toString(this.depths[j])+
								"_p"+ Integer.toString(p) +"\n");
						
						String tmp =" -x 0.5 -t 1 ";
						if (this.depths[j] == 250) {
							tmp =" -x 1.25 -t 2 ";
						}
						
						bw2.write("perl /gpfs/home/chencao.new/LRSIM/simulateLinkedReads.pl -r  $prefix_fa\\.fa  -p"
								+ "  $prefix\\/freq_"+ Integer.toString(this.freq_cutoff[f])
								+ "_pool_"+ Integer.toString(this.pools[i])+ 
								"_dep_"+ Integer.toString(this.depths[j])+
								"_p"+ Integer.toString(p) 
								+ " -c /gpfs/home/chencao.new/LRSIM/test/fragmentSizesList -e 0.001 -E 0.001 "
										+ "  -o -1 90000000  -4 0 -7 0 -0 0  -f 50  -m 10 "+ tmp+ "\n");
						
						
						
						
					}
					bw2.close();

//-------------------------------------------------------------------------------------------------------------------------	
					
					BufferedWriter bw3 = new BufferedWriter(new FileWriter(this.prefix_folder + "/cmd/freq_"+
							Integer.toString(this.freq_cutoff[f])+"_pool_"+  
							Integer.toString(this.pools[i])
							+ "_dep_"+ Integer.toString(this.depths[j]) +".3.cmd"));
					
					bw3.write("#!/bin/bash\n");
					bw3.write("#SBATCH --job-name=f_"+Integer.toString(this.freq_cutoff[f])+
							"_p_"+ Integer.toString(this.pools[i])
							+ "_d_"+ Integer.toString(this.depths[j]) +"\n");
					bw3.write("#SBATCH --workdir="+ this.prefix_folder + "/freq_"+
							Integer.toString(this.freq_cutoff[f])+ "_pool_"+ Integer.toString(this.pools[i])
					+ "_dep_"+ Integer.toString(this.depths[j]) +"\n");
					bw3.write("#SBATCH --error="+this.prefix_folder+ "/cmd/f_"+ Integer.toString(this.freq_cutoff[f])+
							"_p_"+Integer.toString(this.pools[i])
					+ "_d_"+ Integer.toString(this.depths[j])+".3.error\n" );
					bw3.write("#SBATCH --output="+this.prefix_folder+ "/cmd/f_"+ Integer.toString(this.freq_cutoff[f])
							+ "_p_"+Integer.toString(this.pools[i])
					+ "_d_"+ Integer.toString(this.depths[j])+".3.output\n" );
					
					bw3.write("#SBATCH --ntasks=1\n");
					bw3.write("#SBATCH --cpus-per-task=12\n");
					bw3.write("#SBATCH --time=99-00:00:00\n");
					bw3.write("#SBATCH --nodes=1\n");
					
					bw3.write("java=/export/home/jhe/download/java_jdk_8u201/jdk1.8.0_201/bin/java\n");
				
					bw3.write("ref="+ this.genome_path+"\n");
					bw3.write("gatk=/export/qlong/chencao/Work/poolhapx/software/gatk-4.0.0.0/"
							+ "gatk-package-4.0.0.0-local.jar\n");
					
					bw3.write("java=/export/home/jhe/download/java_jdk_8u201/jdk1.8.0_201/bin/java\n");
					
//					/export/qlong/PoolHapX/fix_fq.py
					
//					/export/qlong/PoolHapX/longranger-2.2.2/longranger
					
					
					for (int p=0;p < this.pools[i];p ++) {
						
						bw3.write("prefix="+ this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])+ 
								"_pool_"+ Integer.toString(this.pools[i])+ 
								"_dep_"+ Integer.toString(this.depths[j])+"/input/fastq/freq_"+ Integer.toString(this.freq_cutoff[f])
								+ "_pool_"+ Integer.toString(this.pools[i])+ 
								"_dep_"+ Integer.toString(this.depths[j])+
								"_p"+ Integer.toString(p) +"\n");
						
						bw3.write("cd  $prefix\n");
						
						bw3.write("python  /export/qlong/PoolHapX/fix_fq.py  $prefix\n");
						
						
						bw3.write("/export/qlong/PoolHapX/longranger-2.2.2/longranger  align --id=freq_" + 
								Integer.toString(this.freq_cutoff[f])
						+ "_pool_"+ Integer.toString(this.pools[i])+ 
						"_dep_"+ Integer.toString(this.depths[j])+
						"_p"+ Integer.toString(p)+ 
						"  --reference=/export/qlong/chencao/Work/poolhapx/slim/ref/plasmodium/refdata-plasmodium "
						+ " --fastqs=$prefix\n");
						
						
						
						bw3.write("prefix="+ this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])+ 
								"_pool_"+ Integer.toString(this.pools[i])+ 
								"_dep_"+ Integer.toString(this.depths[j])+"/input/fastq/freq_"+ Integer.toString(this.freq_cutoff[f])
								+ "_pool_"+ Integer.toString(this.pools[i])+ 
								"_dep_"+ Integer.toString(this.depths[j])+
								"_p"+ Integer.toString(p) +"\n");
						
						bw3.write("prefix_bam="+ this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])
								+ "_pool_"+ Integer.toString(this.pools[i])+ 
									"_dep_"+ Integer.toString(this.depths[j])+"/input/fastq/freq_" +Integer.toString(this.freq_cutoff[f])
									+ "_pool_"+ Integer.toString(this.pools[i])+ 
									"_dep_"+ Integer.toString(this.depths[j])+
									"_p"+ Integer.toString(p) +"/freq_"+ Integer.toString(this.freq_cutoff[f])
									+ "_pool_"+ Integer.toString(this.pools[i])+ 
									"_dep_"+ Integer.toString(this.depths[j])+
									"_p"+ Integer.toString(p) +
									"\n");

													
//						bw3.write("gunzip   $prefix\\.bwa.read1.fastq\n");
//						bw3.write("gunzip   $prefix\\.bwa.read2.fastq\n");
//						bw3.write("$bwa mem -t 4 $ref $prefix\\.bwa.read1.fastq $prefix\\.bwa.read2.fastq "
//								+ "| samtools view -Shub - > $prefix_bam\\.bam\n");
//						bw3.write("samtools sort -o  " +	"$prefix_bam\\.srt.bam  $prefix_bam\\.bam\n");
						
						
						bw3.write("inbam="+ this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])
								+ "_pool_"+ Integer.toString(this.pools[i])+ 
										"_dep_"+ Integer.toString(this.depths[j])+"/input/bam/freq_" + Integer.toString(this.freq_cutoff[f])
										+ "_pool_"+ Integer.toString(this.pools[i])+ 
										"_dep_"+ Integer.toString(this.depths[j])+
										"_p"+ Integer.toString(p) +".rg.bam\n");
									
						bw3.write("$java -jar $gatk  AddOrReplaceReadGroups -I  $prefix_bam\\/outs/possorted_bam.bam -O $inbam"
								+ " -R $ref -ID " +  "freq_"+Integer.toString(this.freq_cutoff[f]) +  "_pool_"+ Integer.toString(this.pools[i])+ 
								"_dep_"+ Integer.toString(this.depths[j])+ "_p"+ Integer.toString(p)
								+" -LB NPD -PL Illumina -PU NPD -SM freq_"+ Integer.toString(this.freq_cutoff[f])
								+"_pool_"+ Integer.toString(this.pools[i])+ 
								"_dep_"+ Integer.toString(this.depths[j])+ "_p"+ Integer.toString(p)+ "\n");
						
						bw3.write("samtools index $inbam\n");
											
						bw3.write("prefix_vcf="+ this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])
								+ "_pool_"+ Integer.toString(this.pools[i])+ 
								"_dep_"+ Integer.toString(this.depths[j])+"/input/vcf/freq_"+ Integer.toString(this.freq_cutoff[f])
								+ "_pool_"+ Integer.toString(this.pools[i])+ 
								"_dep_"+ Integer.toString(this.depths[j])+
								"_p"+ Integer.toString(p) +"\n");
						
						
						bw3.write("outgvcf=$prefix_vcf\\.raw.g.vcf\n");
						bw3.write("$java -jar $gatk  HaplotypeCaller -R  $ref -I $inbam"
								+ " -ERC  GVCF -ploidy 8 --heterozygosity 0.01  --max-alternate-alleles 1 -O "
								+ "$outgvcf\n");
					}
	//	Step 4: Join all pool-specific gVCFs into a joint gVCF file and convert to VCF.
					
					bw3.write("prefix="+ this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])
							+ "_pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j])+"/input/vcf/freq_"+Integer.toString(this.freq_cutoff[f])
							+ "_pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j]) +"\n");
					
					String tmp= "$java -jar $gatk CombineGVCFs -R $ref ";
					for (int p=0;p < this.pools[i];p ++) {
						tmp=tmp+" -V " + "$prefix\\"+
								"_p"+ Integer.toString(p)+".raw.g.vcf ";
					}
					
					
					tmp=tmp+ " -O   $prefix\\.g.vcf\n"  ;
					bw3.write(tmp );
					
					bw3.write("$java -jar $gatk   GenotypeGVCFs -R $ref -V  $prefix\\.g.vcf" + 
							" -ploidy 8 -O $prefix\\.raw.vcf\n");
					
					bw3.write("prefix_vcf="+ this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])
							+ "_pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j])+"/input/freq_"+Integer.toString(this.freq_cutoff[f])
							+ "_pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j]) +"\n");
					
					bw3.write( "$java -jar $gatk  SelectVariants -R $ref -V  $prefix\\.raw.vcf" + 
							" -O $prefix_vcf\\.vcf\n");
					
	//	Step 5. Convert each pool-specific BAM file to SAM (i.e.: text), then VEF files.	
					for (int p=0;p < this.pools[i];p ++) {
						
						bw3.write("prefix_bam="+ this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])
						+ "_pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j])+"/input/fastq/freq_" +Integer.toString(this.freq_cutoff[f])
							+ "_pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j])+
							"_p"+ Integer.toString(p) +"/freq_"+ Integer.toString(this.freq_cutoff[f])
							+ "_pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j])+
							"_p"+ Integer.toString(p) +
							"\n");
						
						
						bw3.write("prefix_sam="+this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])
								+ "_pool_"+ Integer.toString(this.pools[i])+ 
								"_dep_"+ Integer.toString(this.depths[j])+"/input/sam/freq_"+Integer.toString(this.freq_cutoff[f])
								+ "_pool_"+ Integer.toString(this.pools[i])+ 
								"_dep_"+ Integer.toString(this.depths[j])+
								"_p"+ Integer.toString(p) +"\n");
						
						
						
						bw3.write("samtools view -ho $prefix_sam\\.sam  $prefix_bam\\/outs/possorted_bam.bam\n");
						
					}
					bw3.write("properties="+this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])
						+ "_pool_"+ Integer.toString(this.pools[i])+ 
						"_dep_"+ Integer.toString(this.depths[j])+"/input/PHX.properties\n");
					bw3.write("poolhapx=/export/qlong/PoolHapX/PoolHapX.jar\n");
					bw3.write("$java -jar $poolhapx format $properties\n");
					bw3.write("$java -jar $poolhapx gc $properties\n");
					bw3.write("$java -jar $poolhapx aem $properties\n");
					bw3.write("$java -jar $poolhapx evaluate $properties\n");
					bw3.write("\n");
					bw3.close();
					
//-------------------------------------------------------------------------------------------------------------------------	
					
					
					BufferedWriter bw4 = new BufferedWriter(new FileWriter(this.prefix_folder + "/cmd/freq_"+
							Integer.toString(this.freq_cutoff[f])+"_pool_"+  
							Integer.toString(this.pools[i])
							+ "_dep_"+ Integer.toString(this.depths[j]) +".4.cmd"));
					
					bw4.write("#!/bin/bash\n");
					bw4.write("#SBATCH --job-name=f_"+Integer.toString(this.freq_cutoff[f])+
							"_p_"+ Integer.toString(this.pools[i])
							+ "_d_"+ Integer.toString(this.depths[j]) +"\n");
					bw4.write("#SBATCH --workdir="+ this.prefix_folder + "/freq_"+
							Integer.toString(this.freq_cutoff[f])+ "_pool_"+ Integer.toString(this.pools[i])
					+ "_dep_"+ Integer.toString(this.depths[j]) +"\n");
					bw4.write("#SBATCH --error="+this.prefix_folder+ "/cmd/f_"+ Integer.toString(this.freq_cutoff[f])+
							"_p_"+Integer.toString(this.pools[i])
					+ "_d_"+ Integer.toString(this.depths[j])+".4.error\n" );
					bw4.write("#SBATCH --output="+this.prefix_folder+ "/cmd/f_"+ Integer.toString(this.freq_cutoff[f])
							+ "_p_"+Integer.toString(this.pools[i])
					+ "_d_"+ Integer.toString(this.depths[j])+".4.output\n" );
					
					bw4.write("#SBATCH --ntasks=1\n");
					bw4.write("#SBATCH --cpus-per-task=12\n");
					bw4.write("#SBATCH --time=99-00:00:00\n");
					bw4.write("#SBATCH --nodes=1\n");
					
					bw4.write("java=/export/home/jhe/download/java_jdk_8u201/jdk1.8.0_201/bin/java\n");
				
					bw4.write("ref="+ this.genome_path+"\n");
					bw4.write("gatk=/export/qlong/chencao/Work/poolhapx/software/gatk-4.0.0.0/"
							+ "gatk-package-4.0.0.0-local.jar\n");
					
					bw4.write("java=/export/home/jhe/download/java_jdk_8u201/jdk1.8.0_201/bin/java\n");
					
//					/export/qlong/PoolHapX/fix_fq.py
					
//					/export/qlong/PoolHapX/longranger-2.2.2/longranger
					
//	Step 5. Convert each pool-specific BAM file to SAM (i.e.: text), then VEF files.	
					for (int p=0;p < this.pools[i];p ++) {
						bw4.write("prefix_bam="+ this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])
						+ "_pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j])+"/input/fastq/freq_" +Integer.toString(this.freq_cutoff[f])
							+ "_pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j])+
							"_p"+ Integer.toString(p) +"/freq_"+ Integer.toString(this.freq_cutoff[f])
							+ "_pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j])+
							"_p"+ Integer.toString(p) +
							"\n");
						
						
						bw4.write("prefix_sam="+this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])
								+ "_pool_"+ Integer.toString(this.pools[i])+ 
								"_dep_"+ Integer.toString(this.depths[j])+"/input/sam/freq_"+Integer.toString(this.freq_cutoff[f])
								+ "_pool_"+ Integer.toString(this.pools[i])+ 
								"_dep_"+ Integer.toString(this.depths[j])+
								"_p"+ Integer.toString(p) +"\n");
						
						
						
						bw4.write("samtools view -ho $prefix_sam\\.sam  $prefix_bam\\/outs/possorted_bam.bam\n");
						
					}
					bw4.write("properties="+this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])
						+ "_pool_"+ Integer.toString(this.pools[i])+ 
						"_dep_"+ Integer.toString(this.depths[j])+"/input/PHX.properties\n");
					bw4.write("poolhapx=/export/qlong/PoolHapX/PoolHapX.jar\n");
					
					bw4.write("$java -jar $poolhapx format $properties\n");
					bw4.write("$java -jar $poolhapx gc $properties\n");
					bw4.write("$java -jar $poolhapx aem $properties\n");
					bw4.write("$java -jar $poolhapx evaluate $properties\n");
					bw4.write("\n");
					bw4.close();
					
					
        	
					
					new File(this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])
							+ "_pool_"+ Integer.toString(this.pools[i])+ 
								"_dep_"+ Integer.toString(this.depths[j])+"/input/").mkdir();         
					
					BufferedWriter bw_properties = new BufferedWriter(new FileWriter(
						this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])
								+ "_pool_"+ Integer.toString(this.pools[i])+ 
						"_dep_"+ Integer.toString(this.depths[j])+"/input/PHX.properties"));
					
					bw_properties.write("Proj_Name = "+ "freq_"+Integer.toString(this.freq_cutoff[f])
							+ "_pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j])+"\n" );
	
					bw_properties.write("Input_Dir = "+this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])
							+ "_pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j]) +"/input/\n" );
					
					bw_properties.write("Intermediate_Dir = "+this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])
							+ "_pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j]) +"/intermediate/\n" );
					
					bw_properties.write("Output_Dir = "+this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])
							+ "_pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j]) +"/output/\n" );
					
					bw_properties.write("Gold_Dir = "+this.prefix_folder + "/freq_"+Integer.toString(this.freq_cutoff[f])
							+ "_pool_"+ Integer.toString(this.pools[i])+ 
							"_dep_"+ Integer.toString(this.depths[j]) +"/gold_standard/\n" );
					
					bw_properties.write("Num_Pos_Window = 1000\n");
					bw_properties.write("Num_Gap_Window = 2\n");
					bw_properties.write("Num_Pos_Job = 1000\n");
					bw_properties.write("In-pool_Gap_Support_Min = 1\n");
					bw_properties.write("All-pool_Gap_Support_Min = 1\n");
					bw_properties.write("Level_1_Region_Size_Min = 11\n");
					bw_properties.write("Level_1_Region_Size_Max = 13\n");
					bw_properties.write("Level_2_Region_Size_Min = 11\n");
					bw_properties.write("Level_2_Region_Size_Max = 13\n");
					bw_properties.write("Est_Ind_PerPool = 1000000\n");
					bw_properties.write("Level_3_4_Region_Mismatch_Tolerance = 1\n");
					bw_properties.write("Level_5_6_Region_Mismatch_Tolerance = 2\n");
					bw_properties.write("Level_7_8_Region_Mismatch_Tolerance = 5\n");
					bw_properties.write("AEM_Maximum_Level = 7\n");
					bw_properties.write("BFS_Mismatch_Tolerance = 8\n");
					bw_properties.write("AEM_Iterations_Max = 200\n");
					bw_properties.write("AEM_Convergence_Cutoff = 0.0\n");
					bw_properties.write("AEM_Zero_Cutoff = 0.0\n");
					bw_properties.write("AEM_Regional_Cross_Pool_Freq_Cutoff = 0.0\n");
					bw_properties.write("AEM_Regional_HapSetSize_Max = 50\n");
					bw_properties.write("AEM_Regional_HapSetSize_Min = 3\n");
					bw_properties.write("Virtual_Coverage_Link_GC = 1500\n");
					bw_properties.write("Hc_Similarity_Cutoff =0.95\n");
					bw_properties.write("MCC_Freq_Cutoff = 0.01\n");
					bw_properties.write("Rscript_path = /export/home/chencao/miniconda2/bin/Rscript\n");
					bw_properties.write("Regression_Distance_Max_Weight = 2.5\n");
					bw_properties.write("Regression_Coverage_Weight = 1.0\n");
					bw_properties.write("Regression_Mismatch_Tolerance= 7\n");
					bw_properties.write("Regression_One_Vector_Weight = 5.0 \n");
					bw_properties.write("Regression_Hap_VC_Weight = 2.0 \n");
					bw_properties.write("Regression_Hap_11_Weight = 1.0 \n");
					bw_properties.write("Regression_Regional_HapSetSize_Max = 15\n");
					bw_properties.write("Regression_Regional_HapSetSize_Min = 40\n");
					bw_properties.write("Regression_Gamma_Min = 0.0001\n");
					bw_properties.write("Regression_n_Gamma = 10\n");
					bw_properties.write("Regression_Gamma_Max = 0.1\n");
					bw_properties.write("Regression_Maximum_Regions = 2\n");
					bw_properties.write("Maximum_Selected_HapSet = 25\n");
					bw_properties.write("Sequencing_Technology = 10X_linked_reads\n");
					bw_properties.write("Number_Threads = 3\n");		
					bw_properties.write("Species = 10X\n");	 
					bw_properties.close();
					
				}
			}
		}
		
	}

	public static  void main(String[] args) throws IOException, InterruptedException {
		System.out.println("PoolHapX Comparison Simulation... ...");
		ComparisonSim_10X cs = new ComparisonSim_10X(args[0]);
		System.out.println("Done, Enjoy!");
	}	
}


