package MiscFunctions;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import PoolHap.Entrance;

public class BAMFormatterGATK {
	
	/**
	 * @author quanlong 2019-07
	 * 
	 * Generate vef files and var_freq file based on SAM files and a single VCF file that contains 
	 * variants information. 
	 * 
	 * Before entering this function, Entrance.get_filepaths() will sort out the Pool-IDs and their
	 * orders by reading the SAM files in the input_dir/sam/ folder. This order will be fixed 
	 * throughout the PoolHapX program.   
	 * 
	 * @param input_dir
	 * @param inter_dir
	 * @param project_name
	 * @param pool_IDs
	 * @throws IOException
	 */
	
	String sequencing_technology;
	
	public static void rewrite_inter_vars(String gold_stand_vars, String inter_vars) throws IOException {
		
		BufferedReader VarsReader = new BufferedReader(new FileReader(inter_vars));
		String currLine ="";
		String header="";
		HashMap<String , String> pos_dict = new HashMap<String, String>(); 
		int count =0;
		while ((currLine = VarsReader.readLine()) != null) {
			currLine = currLine.replace(" ", "");
			if (count ==0) {
				header =currLine;
			}else {
				String tmp = currLine;
				String[] fullLine = currLine.split("\t");
				pos_dict.put(fullLine[0], tmp); 
			}
			count++;
			
		}
		VarsReader.close();
		count =0;
//		for (String entry : pos_dict.keySet()) {
//			System.out.println(pos_dict.get(entry));
//		}
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(inter_vars));
		bw.write(header+"\n");
		BufferedReader GoldReader = new BufferedReader(new FileReader(gold_stand_vars));
		
		while ((currLine = GoldReader.readLine()) != null) {
			currLine = currLine.replace(" ", "");
			if (count>0) {
				String[] fullLine = currLine.split("\t");
//				System.out.println(fullLine[0]);
				if (pos_dict.containsKey(fullLine[0])) {
					bw.write(pos_dict.get(fullLine[0])+"\n");
				} else {
					String tmp = fullLine[0];
					for (int i=1; i< fullLine.length; i++) {
						tmp =tmp+"\t"+  "0";
					}
					bw.write(tmp+"\n");
				}
			}
			count++;
		}
		
		GoldReader.close();
		bw.close();
		
		return ;
	}
			
	
	public static void generate_vef_varfreq(String input_dir, String inter_dir, 
			String project_name, String[] sam_files, 
			HashMap<String, Integer> pool_name2index, String seq_tech) throws IOException {
		
		Boolean simmode = false; // Boolean.parseBoolean(args[3]);
		HashSet<Integer> trueVarPos = new HashSet<Integer>(); 
		int numPts=pool_name2index.size();
		String[] pool_IDs=new String[numPts];
		for(String pool_id: pool_name2index.keySet()) {
			int index=pool_name2index.get(pool_id);
			pool_IDs[index]=pool_id;
		}
		if (simmode) {
			trueVarPos = simvarsReader(inter_dir + "/simvars.intra_freq.txt", numPts); 
		}
		//PrintWriter VEFList = new PrintWriter(inter_dir + project_name + ".vef.list"); 
		HashMap<Integer,HashMap<String,VarObj>> variantEncoder = VariantMapping(
				input_dir + project_name + ".vcf", true, 
				inter_dir + project_name+"_vars.intra_freq.txt", 
				pool_name2index, simmode, trueVarPos);	
		// args[1] needs to be interdir + prefix + .all.vcf.
		new File(inter_dir+"/vef_raw/").mkdir();
		new File(inter_dir+"/vef/").mkdir();
		for (int p = 0; p < numPts; p++) {
			String raw_vef=inter_dir+"/vef_raw/" + pool_IDs[p] + ".raw.vef";
			String linked_vef=inter_dir+"/vef/" + pool_IDs[p] + ".vef";
//			System.out.println(seq_tech);
			VEFMaker(sam_files[p], raw_vef, variantEncoder, seq_tech);
			PairedReadLinker.link_paired_vef(raw_vef, linked_vef);
		//	VEFList.println(outdir + prefix + "_p" + p + ".vef"); 
		}
//		VEFList.close();
	}
	
	public static HashSet<Integer> simvarsReader(String simvarsFile, int numPts) throws IOException {
		BufferedReader SVReader = new BufferedReader(new FileReader(simvarsFile)); 
		String currLine = SVReader.readLine();	// Skip the first line.
		currLine = SVReader.readLine();
		HashSet<Integer> trueVarPos = new HashSet<Integer>();  
		svRead: while (currLine != null) {
			String[] fullLine = currLine.split("\t"); // 0;211;211;0:1		0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
			for (int p = 2; p <= numPts + 1; p++) {
				if (Double.parseDouble(fullLine[p]) == 0)  {
					currLine = SVReader.readLine();
					continue svRead;
				}
			}
			int truePos = Integer.parseInt(fullLine[0].split(";")[1]); 
			trueVarPos.add(truePos); 	
			currLine = SVReader.readLine();
		}
		SVReader.close();
		return trueVarPos;
	}

	public static HashMap<Integer,HashMap<String,VarObj>> VariantMapping(
			String input_VCF, boolean biallelic, 
			String out_put_vars_freq_file, 
			HashMap<String, Integer> pool_name2index, Boolean simmode, 
			HashSet<Integer> trueVarPos) throws IOException { 	
		int numPts=pool_name2index.size();
		// name2vcf_index maps the names in SAM files to the columns in the multi-individual VCF file. 
		HashMap<String, Integer> name2vcf_index=new HashMap<String, Integer>();
		BufferedReader VCFReader = new BufferedReader(new FileReader(input_VCF));
		String currLine = VCFReader.readLine();
		// System.out.println(currLine);
		while (!currLine.contains("#CHROM")) currLine = VCFReader.readLine();
		String[] header_line = currLine.split("\t"); // Match up the order of the pool IDs in the VCF to their actual row number.
		int[] actual_row = new int[numPts]; 
		for (int i = 9; i < header_line.length; i++) {
			//actual_row[i - 9] = Integer.parseInt(header_line[i]); // Lauren's old code
			actual_row[i - 9] = pool_name2index.get(header_line[i]); 
			// above uses new mapping using actual pool-IDs from SAM files 
		}
		currLine = VCFReader.readLine();
		// variantEncoder stores position to alleles mapping.
		HashMap<Integer,HashMap<String,VarObj>> variantEncoder = new HashMap<Integer,HashMap<String,VarObj>>();
		ArrayList<double[]> poolFreqs = new ArrayList<double[]>();
		ArrayList<Integer> posTracker = new ArrayList<Integer>();
		int variant_index = 0; 
		while (currLine != null) {
			String[] var_info = currLine.split("\t");
			int pos = Integer.parseInt(var_info[1]);
			if (simmode)
				if (!trueVarPos.contains(pos)) {
					currLine = VCFReader.readLine();
					continue;	// If this is not a true variant position, skip it.
				}
			posTracker.add(pos);
			variantEncoder.put(pos, new HashMap<String,VarObj>());
			// var_info[4] is the 5th column in a VCF file, denoting the ALT (alternative alleles).
			Scanner altAlleleScanner = new Scanner(var_info[4]);	
			// This Scanner runs through the String containing all of the alternate allele(s). There is at least one. 
			altAlleleScanner.useDelimiter(",");
			int variantCode = 1;
			while(altAlleleScanner.hasNext()) {					// In case of highly polymorphic (>=2 alternate alleles) positions.
				String variantNow = altAlleleScanner.next().replace("\"","");
				VarObj altAlleleAtPos = new VarObj(pos,variantCode);	// This is the information that will be reported when the alternate allele at this position in the read is accessed.
				variantEncoder.get(pos).put(variantNow,altAlleleAtPos);			// The information is mapped to the alternate allele. 
				variantCode++;										// The code corresponding to the alternate allele is incremented so no two alternate alleles in the same location have the same code.
				if (biallelic == true) break; 	// If PoolHap has not been adjusted for 3+ alleles. 
			} altAlleleScanner.close();
			poolFreqs.add(new double[numPts]); 
			for (int p = 0; p < numPts; p++) {
				String[] varCts = var_info[p + 9].split(":")[1].split(",");	
				// This is the AD of the GT:AD:DP:GQ:PL of the p0 block of genotypes. 
				double ref = Double.parseDouble(varCts[0]); 
				double alt = Double.parseDouble(varCts[1]); 
				if (ref != 0 || alt != 0) poolFreqs.get(variant_index)[actual_row[p]] = alt / (ref + alt);
				else poolFreqs.get(variant_index)[actual_row[p]] = 0;	// In case there are no reads at this position at all.
				// System.out.println(poolFreqs.get(v)[p]);
			}
			variant_index++; 
			currLine = VCFReader.readLine();
		}
		VCFReader.close();
		System.out.println("Reading VCF file: Done.");
		
		// Write the in-pool variants frequencies to file. The IDs of the pools are based on the 
		// pre-specified order in the pool_name2index;
		String[] pool_IDs=new String[numPts];
		for(String pool_id: pool_name2index.keySet()) {
			int index=pool_name2index.get(pool_id);
			pool_IDs[index]=pool_id;
		}
		BufferedWriter br = new BufferedWriter(new FileWriter(out_put_vars_freq_file));
		br.write("Pool_ID"); 
		for(int p=0;p<numPts;p++) br.write("\t"+pool_IDs[p] );
		for(int a=0;a<variant_index;a++) {
			br.write("\n0;" + posTracker.get(a) + ";" + posTracker.get(a) + ";0:1"); 
			// TODO This only allows for biallelic simple loci (single alternate allele) for now. 
			for(int p=0;p<numPts;p++) br.write("\t" + poolFreqs.get(a)[p]);
		}
		br.close();
		System.out.println("In-Pool variants frequncy have been written to "+ out_put_vars_freq_file);
		return variantEncoder;
	}
	
	public static void VEFMaker(String input_SAM, String output_vef_raw,
			HashMap<Integer,HashMap<String,VarObj>> variantEncoder, String seq_tech) 
			throws FileNotFoundException {
		Set<Integer> keyMATCH, keyINS, keyDEL, variantKS = variantEncoder.keySet(), tempSet, indelKS; 
		SortedSet<Integer> variantSS = new TreeSet<Integer>();
		variantSS.addAll(variantKS); 

		File BAMPath = new File(input_SAM); 
		Scanner BAMScanner = new Scanner(BAMPath);
		BAMScanner.useDelimiter("\t");
		String currLine = BAMScanner.nextLine(); 
		// System.out.println(currLine);
		while (currLine.matches("@(.*)")) {
			// String prevLine = currLine;
			currLine = BAMScanner.nextLine(); 
			// if (prevLine.contains("PN:")) break;
			// System.out.println(currLine);
		}
		String QNAME, CIGAR, SEQ, tempPos = "",currAltAllele;
		Integer FLAG, startPOS = 1, currentPos, posToAdd, matchToIndel, skipSBases = 0;
		VarObj currVarObj; 
		HashMap<Integer, HashMap<Integer,String>> hmMATCH = new HashMap<Integer, HashMap<Integer,String>>();
		HashMap<Integer, HashMap<Integer,String>> hmINS = new HashMap<Integer, HashMap<Integer,String>>();
		HashMap<Integer, HashMap<Integer,String>> hmDEL = new HashMap<Integer, HashMap<Integer,String>>();
		HashMap<Integer,String> defaultHM, tempFinderHM;
		HashMap<String,VarObj> tempReporterHM; 
		/* HashMap<Integer,Integer> hmVarCount = new HashMap<Integer,Integer>(); 
		for (Integer i : variantSS) {
			hmVarCount.put(i,0);
		}
		int readCount = 0; */
		//PrintWriter VEFFile = new PrintWriter(BAMPrefix + ".raw.vef");// replaced with the line below
		PrintWriter VEFFile = new PrintWriter(output_vef_raw);
		while (BAMScanner.hasNextLine()) {
			// readCount++; 
			QNAME = BAMScanner.next();
			String QNAME_tmp = QNAME;
			if (seq_tech.equals("10x_linked-reads")) {
//				System.out.println(seq_tech);
				QNAME="";
			}
			// System.out.println(QNAME);
			FLAG = BAMScanner.nextInt();
			if (FLAG > 163) {
				BAMScanner.nextLine();
				continue; 	// If the read is chimeric, or mapped to more than one location, or of poor technical quality, skip it. 
			}
			BAMScanner.next(); 	// Skip the FLAG, RNAME.
			startPOS = BAMScanner.nextInt();
			currentPos = startPOS; 
			BAMScanner.next(); 		// Skip the MAPQ.			
			CIGAR = BAMScanner.next();
			
			if (CIGAR.charAt(0) == '*') {
				BAMScanner.nextLine(); 	// * refers to the fact that no CIGAR information is available. 
				continue;	 			// The above line skips the rest of the information in this read entirely.
			}
			for (int c = 0; c < CIGAR.length(); c++) {
				Character tempChar = CIGAR.charAt(c); 
				if (tempChar.compareTo('M') == 0) {
					posToAdd = Integer.parseInt(tempPos);
					for (int addPos = currentPos; addPos < currentPos + posToAdd; addPos++) {
						defaultHM = new HashMap<Integer,String>(); 
						defaultHM.put(1, "N"); 
						hmMATCH.put(addPos,defaultHM);
					}
					currentPos += posToAdd;
					tempPos = ""; 
				} else if (tempChar.compareTo('I') == 0) { 
					// Command to check for insertions: samtools view HIV_p1_sample_1.procd.bam | awk '{if ($6 ~ /I/) print $6;}' -
					posToAdd = Integer.parseInt(tempPos);
					matchToIndel = currentPos - 1; 	// This is the first position of the insertion, which it is recognized by. 
					hmMATCH.remove(matchToIndel); 	// Remove the first position of the insertion from the match HM.
					defaultHM = new HashMap<Integer,String>(); 
					defaultHM.put(posToAdd + 1, "N");	// This accounts for the entire length of the insertion.
					hmINS.put(matchToIndel, defaultHM); 
					currentPos = matchToIndel + 1;		// Since the insertion starts at the previous currentPos - 1, the next match along starts at matchToIndel + 1 = previous currentPos.
					tempPos = ""; 
				} else if (tempChar.compareTo('D') == 0) {
					posToAdd = Integer.parseInt(tempPos);
					matchToIndel = currentPos - 1; 	// This is the first position of the deletion, which it is recognized by. 
					hmMATCH.remove(matchToIndel); 	// Remove the first position of the deletion from the match HM.
					defaultHM = new HashMap<Integer,String>(); 
					defaultHM.put(posToAdd + 1, "N");	// This accounts for the entire length of the deletion.
					hmDEL.put(matchToIndel, defaultHM); 
					currentPos += posToAdd; 	// Since the deletion starts at the previous currentPos - 1, the next match along starts at currentPos + matchToIndel - 1.
					tempPos = ""; 
				} else if (tempChar.compareTo('S') == 0) {
					skipSBases = Integer.parseInt(tempPos);	// These bases do appear in SEQ BUT they are not aligned to the reference genome. Therefore these bases do not contribute to the currentPOS. 
					tempPos = ""; 				
				} else if (tempChar.compareTo('H') == 0) {
					tempPos = ""; // These bases DO NOT appear in SEQ and they are not aligned to the reference genome. Therefore these bases do not contribute to the currentPOS. 				
				} else { 
					tempPos += tempChar;
				}
			}

			keyMATCH = hmMATCH.keySet(); 
			keyINS = hmINS.keySet();
			keyDEL = hmDEL.keySet(); 
			for (int s = 0; s < 3; s++) {
				BAMScanner.next(); 	// Skip the RNEXT, PNEXT, TLEN.
			}
			SEQ = BAMScanner.next();
			int endPOS = startPOS + SEQ.length() - 1;
			currentPos = startPOS; 
			// 4. Test the following to see if the SEQ-parsing code works properly, and the proper bases are assigned to the correct Hash Maps. 
			for (int s = skipSBases; s < SEQ.length(); s++) {	// Skip all of the bases that were soft-clipped but still reported at the beginning of SEQ. In the event there is no 'S' in the CIGAR string, this is 0 and we start reporting from the end of SEQ.  
				// System.out.println(s + "\t" + currentPos);
				if (keyMATCH.contains(currentPos)) {
					hmMATCH.get(currentPos).put(1,SEQ.substring(s,s+1));
					currentPos++; 
				} else if (keyINS.contains(currentPos)) {
					tempSet = hmINS.get(currentPos).keySet();
					for (Integer length : tempSet) {
						hmINS.get(currentPos).put(length,SEQ.substring(s,s+length)); 
					}
					currentPos++; 
				} else if (keyDEL.contains(currentPos)) {
					tempSet = hmDEL.get(currentPos).keySet();
					for (Integer length : tempSet) {
						hmDEL.get(currentPos).put(length,SEQ.substring(s,s+1));
						currentPos += length; 
					}
				}
			}
			
			/* for (int i : keyMATCH) {
				System.out.println(i + "\t" + hmMATCH.get(i));
			} */
			skipSBases = 0; // Reset the 'start' position of SEQ-processing to 0 in case there aren't any S-es in the next CIGAR string. 
			
			StringBuilder readInfo = new StringBuilder(); 
			// System.out.println(startPOS + "\t" + endPOS);
			for (Integer VCFVarPos : variantSS.subSet(startPOS,endPOS+1)) {	// subSet method is (inclusive,exclusive)
				// System.out.println(VCFVarPos); 
				tempReporterHM = variantEncoder.get(VCFVarPos);
				if (keyMATCH.contains(VCFVarPos)) {
					currAltAllele = hmMATCH.get(VCFVarPos).get(1);
					// System.out.println(currAltAllele);
					if (tempReporterHM.containsKey(currAltAllele)) {
						currVarObj = tempReporterHM.get(currAltAllele);
						// VEFFile.print(currVarObj.intM + "=" + currVarObj.varCode + ";");
						readInfo.append(currVarObj.intM + "=" + currVarObj.varCode + ";");
						// hmVarCount.put(VCFVarPos, hmVarCount.get(VCFVarPos) + 1);
					} else {
						readInfo.append(VCFVarPos + "=0;"); // Fixed as of 10282018. Noted by CC that in some reads, variant positions in those reads were not being annotated.
						//System.out.println(VCFVarPos + " done");	// Realized that the 'else' clause was in the wrong place i.e.: alleles  =/= alternate (i.e: the reference) 
					}												// were not noted when they existed. 
					// System.out.println("match");
				} else if (keyINS.contains(VCFVarPos)) {
					tempFinderHM = hmINS.get(VCFVarPos);
					indelKS = tempFinderHM.keySet(); 
					for (Integer indel : indelKS) {
						currAltAllele = tempFinderHM.get(indel);
						if (tempReporterHM.containsKey(currAltAllele)) {
							currVarObj = tempReporterHM.get(currAltAllele);
							// VEFFile.print(currVarObj.intM + "=" + currVarObj.varCode + ";");
							readInfo.append(currVarObj.intM + "=" + currVarObj.varCode + ";");
							// hmVarCount.put(VCFVarPos, hmVarCount.get(VCFVarPos) + 1);
						} else {
							readInfo.append(VCFVarPos + "=0;");
							// System.out.println(VCFVarPos + " done");
						}
					}
					// System.out.println("ins");
				} else if (keyDEL.contains(VCFVarPos)) {
					tempFinderHM = hmDEL.get(VCFVarPos);
					indelKS = tempFinderHM.keySet(); 
					for (Integer indel : indelKS) {
						currAltAllele = tempFinderHM.get(indel);
						if (tempReporterHM.containsKey(currAltAllele)) {
							currVarObj = tempReporterHM.get(currAltAllele);
							// VEFFile.print(currVarObj.intM + "=" + currVarObj.varCode + ";");
							readInfo.append(currVarObj.intM + "=" + currVarObj.varCode + ";");
							// hmVarCount.put(VCFVarPos, hmVarCount.get(VCFVarPos) + 1);
						} else {
							readInfo.append(VCFVarPos + "=0;");
							// System.out.println(VCFVarPos + " done");
						}
						// System.out.println("del");
					}
				}
			}
			// System.out.println(readInfo.toString());
			if (readInfo.toString().isEmpty()) {
				BAMScanner.nextLine(); 	// 
				hmMATCH.clear();
				hmINS.clear();
				hmDEL.clear(); 
				continue;
			}
			// System.out.println(readInfo);
			int flag=0;
			if (seq_tech.equals("10x_linked-reads")) {
			

				String DM= "";
				
//				OM:i:60 XM:Z:0  TQ:Z:>@>>=?@    TR:Z:TAGGGTT    AS:i:0  XS:i:-70        XT:i:0  BX:Z:TGTCACCAGGGTACGT-1 
				String tmp ="";
				String BX="";
				
				for (int s = 0; s < 14; s++) {
					tmp = BAMScanner.next(); 	
					if ((tmp.length()>4) && (tmp.substring(0, 3).equals("BX:")) ) {
						BX= tmp;
						flag++;
					}
					if ((tmp.length()>4) && (tmp.substring(0, 3).equals("DM:")) ) {
						DM= tmp;
						flag++;
					}
					if (flag==2) {
						break;
					}
				}
				QNAME=DM.replace("\n", "")+":"+ BX.replace("\n", "");
//				System.out.println(QNAME);
			}
			if (flag==2 ) {
				String[] x_arr =  QNAME_tmp.split("_");
				String hap = x_arr[0]+"_"+ x_arr[1]+ "_"+ x_arr[2];
				VEFFile.println(QNAME+"_"+hap + ":\t" + readInfo +"\t"+QNAME_tmp+ "\t"+ startPOS+"\t" + endPOS);
			} 
			if (!seq_tech.equals("10x_linked-reads")) {
				VEFFile.println(QNAME + ":\t" + readInfo +"\t"+QNAME_tmp+ "\t"+ startPOS+"\t" + endPOS);
			}
			
			BAMScanner.nextLine(); 	// 
			hmMATCH.clear();
			hmINS.clear();
			hmDEL.clear(); 
		}
		VEFFile.close();
		BAMScanner.close();
		System.out.println(output_vef_raw +" has been generated.");
	}
}
class VarObj {

	public Integer intM; 	// The position of the variant according to the metareference.
	public Integer varCode; // The numerical code corresponding to the desired alternative allele.  
	
	public VarObj(Integer i, Integer c) {	// This is the constructor. 
		intM = i;
		varCode = c; 
	}
}