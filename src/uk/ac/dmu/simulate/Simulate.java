package uk.ac.dmu.simulate;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.io.iterator.IteratingSMILESReader;
import org.openscience.cdk.layout.HydrogenPlacer;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;

public class Simulate {

        public static void main(String[] args) throws Exception {
        		Properties props=new Properties();
        		props.load(new FileReader("nmrproc.properties"));
        		String projectdir=new String(props.getProperty("datadir")+File.separator+args[0]+File.separator);
        		File localprops=new File(projectdir+"nmrproc.properties");
        		if(localprops.exists()) {
        			Properties props2=new Properties();
        			props2.load(new FileReader(localprops));
        			for(Object prop : props2.keySet()) {
        				props.setProperty((String)prop, props2.getProperty((String)prop));
        			}
        		}
        		boolean doHsqcTocsy=false;
        		if(props.containsKey("usehsqctocsy") && props.get("usehsqctocsy").equals("true"))
        			doHsqcTocsy=true;
        		boolean doHmbc=true;
        		if(props.containsKey("usehmbc") && props.get("usehmbc").equals("false"))
        			doHmbc=false;
        		boolean dotwobonds=true;
        		if(props.containsKey("dotwobonds") && props.get("dotwobonds").equals("false"))
        			dotwobonds=false;
        		boolean debug=false;
        		if(props.containsKey("debug") && props.get("debug").equals("true"))
        			debug=true;
        		boolean use3d = true;
        		boolean usedeeplearning=false;
        		if(props.containsKey("usedeeplearning") && props.get("usedeeplearning").equals("true"))
        			usedeeplearning=true;
                /* get predictor or readers for respredict prediction*/             
            	PredictionTool predictor = null;
            	JSONArray predsc = null;
            	JSONArray predsh = null;
            	if(!usedeeplearning) {
            		predictor = new PredictionTool();
            	}else {
	                JSONParser parser = new JSONParser();
	                File jsonc=new File(projectdir+"predc.json");
	                FileReader fis=new FileReader(jsonc);
	                Object obj = parser.parse(fis);
	                JSONObject cpred = (JSONObject)obj;
	                predsc=(JSONArray)cpred.get("predictions");
	                parser = new JSONParser();
	                File jsonh=new File(projectdir+"predh.json");
	                fis=new FileReader(jsonh);
	                obj = parser.parse(fis);
	                JSONObject hpred = (JSONObject)obj;
	                predsh=(JSONArray)hpred.get("predictions");
            	}
                FileOutputStream fos=new FileOutputStream(projectdir+File.separatorChar+"result"+File.separatorChar+props.getProperty("predictionoutput"));
                FileOutputStream foshsqc=null;
                FileOutputStream foshmbc=null;
                FileOutputStream foshsqctocsy=null;
                File namesfile=new File(projectdir+props.getProperty("msmsinput").substring(0,props.getProperty("msmsinput").length()-4)+"names.txt");
                if(namesfile.exists()) {
                	foshsqc=new FileOutputStream(new File(projectdir+File.separatorChar+"result"+File.separatorChar+props.getProperty("predictionoutput")+"hsqc"));
                	if(doHmbc)
                		foshmbc=new FileOutputStream(new File(projectdir+File.separatorChar+"result"+File.separatorChar+props.getProperty("predictionoutput")+"hmbc"));
                	if(doHsqcTocsy)
                		foshsqctocsy=new FileOutputStream(new File(projectdir+File.separatorChar+File.separatorChar+"result"+props.getProperty("predictionoutput")+"hsqctocsy"));
                }
                IteratingSMILESReader smilesreader = new IteratingSMILESReader(new FileInputStream(new File(projectdir+props.getProperty("msmsinput"))), DefaultChemObjectBuilder.getInstance());
/* get molecule for reader */
                String solvent=props.getProperty("solvent");
                BufferedReader br=null;
                if(namesfile.exists())
                	br=new BufferedReader(new FileReader(namesfile));
                int k=0;
                while(smilesreader.hasNext()) {
                	JSONArray predsmol=null;
                	JSONArray predsmolh=null;
                	if(usedeeplearning) {
                		predsmol=(JSONArray)((JSONObject)predsc.get(k)).get("shifts");
                		predsmolh=(JSONArray)((JSONObject)predsh.get(k)).get("shifts");
                	}
                	//System.out.println(k++);
                	IAtomContainer mol = smilesreader.next();
                	for(IAtom atom : mol.atoms()) {
                		if(atom.getSymbol().equals("H")) {
                			mol.removeAtomAndConnectedElectronContainers(atom);
                		}
                	}
                    StructureDiagramGenerator sdg = new StructureDiagramGenerator();
                    sdg.setMolecule(mol);
                    sdg.generateCoordinates();
                    mol = sdg.getMolecule();
/* !!! No explicit H in mol !!! */
/* add explicit H atoms */
                	addAndPlaceHydrogens(mol);
                	int a=0;
                	for(IAtom atom : mol.atoms()) {
                		if(atom.getID()==null)
                			atom.setID("a"+a);
                		a++;
                	}
                    //we make a clone and remove non h-bearing atoms to isolate spin systems
                	IAtomContainer ac2=new AtomContainer(mol);
                	for(IAtom atom : ac2.atoms()){
                		if(!atom.getSymbol().equals("H")){
        	        		boolean hashydrogen=false;
        	        		for(IAtom neighbour : ac2.getConnectedAtomsList(atom)){
        	        			if(neighbour.getSymbol().equals("H")){
        	        				hashydrogen=true;
        	        				break;
        	        			}
        	        		}
        	        		if(!hashydrogen)
        	        			ac2.removeAtomAndConnectedElectronContainers(atom);
                		}
                	}
                	IAtomContainerSet spinsystems=ConnectivityChecker.partitionIntoMolecules(ac2);
/* detect aromaticity */
                	CDKHueckelAromaticityDetector.detectAromaticity(mol);
/* current atom and current result */
                	IAtom curAtom;
/* write header */
                	//System.out.println("  C    H    c shift  h shift");
                	//System.out.println(" ---------------------------");
/* loop over atoms in mol */
                	if(namesfile.exists()) {
                		String name=br.readLine()+"\r\n";
                		if(doHmbc)
                			foshmbc.write(name.getBytes());
                		foshsqc.write(name.getBytes());
                		if(doHsqcTocsy)
                			foshsqctocsy.write(name.getBytes());
                		//System.out.println(name);
                	}
                	int predictioncount=0;
                	int spheres=0;
            		//fos.write(new String("\nhsqc\n").getBytes());
                	//we use this map to get atoms by id
                	Map<String, IAtom> atomsbyid = new HashMap<String, IAtom>();
                	List<String> done=new ArrayList<String>();
                	for(int i = 0 ; i < mol.getAtomCount() ; i++){
                        curAtom = mol.getAtom(i);
                        atomsbyid.put(curAtom.getID(), curAtom);
                        if(curAtom.getAtomicNumber() == 6) {
                        	float[] resultc = null;
                        	if(usedeeplearning) {
                        		resultc = shift(predsmol, i);
                        	}else{
                        		resultc = predictor.predict(mol, curAtom, use3d, solvent);
                        	}
                            predictioncount++;
                        	spheres+=resultc[4];
                        	//System.out.println(i+" c "+resultc[4]);
                        	//from a carbon atom, we get the hs which are 2 or 3 bonds away
                        	List<IAtom> away1 = mol.getConnectedAtomsList(curAtom);
                        	for(IAtom away1atom : away1) {
                        		//hsqc
                        		//these are 1 bonds away, if it's a hydrogen, it's a match
                        		if(away1atom.getAtomicNumber()==1  ) {
                                	float[] resulth = null;
                                	if(usedeeplearning) {
                                		resulth = shift(predsmolh, mol.getAtomNumber(away1atom));
                                	}else{
                                		resulth = predictor.predict(mol, away1atom, use3d, solvent);
                                	}
                                	predictioncount++;
                                	spheres+=resulth[4];

                                	//System.out.println(i+" h "+resultc[4]);

                                	if(!done.contains(resultc[1]+";"+resulth[1])) {
                                		//System.out.format(Locale.US, "%3d   %3d%8.2f%8.2f\n", i+1, mol.getAtomNumber(away1atom)+1, resultc[1], resulth[1]);
                                		fos.write(new String(resultc[1]+","+resulth[1]+",q\n").getBytes());
                                		if(debug)
                                			foshsqc.write(new String(resultc[1]+","+resulth[1]+"\n").getBytes());
                                		//System.out.println(resultc[4]+" "+resulth[4]);
                                    	done.add(resultc[1]+";"+resulth[1]);
                                	}
                        		}
                        	}
                        }
                	}
            		//fos.write(new String("hmbc\n").getBytes());
                	done=new ArrayList<String>();
                	List<String> donetocsy=new ArrayList<String>();
                	for(int i = 0 ; i < mol.getAtomCount() ; i++){
                        curAtom = mol.getAtom(i);
                        if(curAtom.getAtomicNumber() == 6) {
                        	//System.out.println(i);
                        	boolean hashydrogen=false;
                        	float[] resultc = null;
                        	if(usedeeplearning) {
                        		resultc = shift(predsmol, i);
                        	}else{
                        		resultc = predictor.predict(mol, curAtom, use3d, solvent);
                        	}
                        	predictioncount++;
                        	spheres+=resultc[4];
                        	//from a carbon atom, we get the hs which are 2 or 3 bonds away
                        	List<IAtom> away1 = mol.getConnectedAtomsList(curAtom);
                        	if(doHmbc) {
	                        	for(IAtom away1atom : away1) {
	                        		if(away1atom.getAtomicNumber()==1)
	                        			hashydrogen=true;
	                            	List<IAtom> away2 = mol.getConnectedAtomsList(away1atom);
	                            	for(IAtom away2atom : away2) {
	                            		if(dotwobonds) {
		                            		//these are 2 bonds away, if it's a hydrogen, it's a match
		                            		if(away2atom.getAtomicNumber()==1 && away1atom.getAtomicNumber() == 6) {
		                                    	float[] resulth = null;
		                                    	if(usedeeplearning) {
		                                    		resulth = shift(predsmolh, mol.getAtomNumber(away2atom));
		                                    	}else{
		                                    		resulth = predictor.predict(mol, away2atom, use3d, solvent);
		                                    	}
		                                    	predictioncount++;
		                                    	spheres+=resulth[4];
		                                    	if(!done.contains(resultc[1]+";"+resulth[1])) {
		                                    		//System.out.format(Locale.US, "%3d   %3d%8.2f%8.2f\n", i+1, mol.getAtomNumber(away2atom)+1, resultc[1], resulth[1]);
		                                    		fos.write(new String(resultc[1]+","+resulth[1]+",b\n").getBytes());
		                                    		if(debug)
		                                    			foshmbc.write(new String(resultc[1]+","+resulth[1]+"\n").getBytes());
		                                    		//System.out.println(resultc[1]+" "+resulth[1]);
		                                        	done.add(resultc[1]+";"+resulth[1]);
		                                    	}
		                            		}
		                            	}
	                                	List<IAtom> away3 = mol.getConnectedAtomsList(away2atom);
	                                	for(IAtom away3atom : away3) {
	                                		//these are 3 bonds away, if it's a hydrogen, it's a match
	                                		if(away3atom.getAtomicNumber()==1 && away2atom.getAtomicNumber() == 6) {
	                                			float[] resulth = null;
	                                			if(usedeeplearning) {
	                                				resulth = shift(predsmolh, mol.getAtomNumber(away3atom));
	                                			}else{
	                                				resulth = predictor.predict(mol, away3atom, use3d, solvent);
	                                			}
		                                    	predictioncount++;
		                                    	spheres+=resulth[4];
		                                    	if(!done.contains(resultc[1]+";"+resulth[1])) {
		                                    		//System.out.format(Locale.US, "%3d   %3d%8.2f%8.2f\n", i+1, mol.getAtomNumber(away3atom)+1, resultc[1], resulth[1]);
		                                    		//System.out.println(resultc[1]+" "+resulth[1]);
		                                    		fos.write(new String(resultc[1]+","+resulth[1]+",b\n").getBytes());
		                                    		if(debug)
		                                    			foshmbc.write(new String(resultc[1]+","+resulth[1]+"\n").getBytes());
		                                        	done.add(resultc[1]+";"+resulth[1]);
		                                    	}
	                                		}
	                                	}
                            		}
                            	}
                        	}
                        	//for hsqctocsy, we look inside spin system
                        	if(doHsqcTocsy && hashydrogen) {
	                        	for(IAtomContainer spinsystem : spinsystems.atomContainers()) {
		                        	for(IAtom atomin2 : spinsystem.atoms()) {
		                        		if(atomin2.getID().equals(curAtom.getID())) {
		                        			//System.out.println("in");
		                                	for(IAtom atomin22 : spinsystem.atoms()) {
		                                		for(IAtom atom3 : mol.atoms()) {
		        	                        		if(atom3.getID().equals(atomin22.getID())) {
		        	                        			//System.out.println("equals");
		        	                        			for(IAtom atom4 : mol.getConnectedAtomsList(atom3)) {
		        	                        				//System.out.println("connected");
		        	                        				if(atom4.getAtomicNumber()==1 && atom3.getAtomicNumber() == 6) {
		        	                        					//System.out.println("1");
		        	                        					float[] resulth = null;
		        	                        					if(usedeeplearning) {
		        	                        						resulth = shift(predsmolh, mol.getAtomNumber(atom4));
		        	                        					}else{
		        	                        						resulth = predictor.predict(mol, atom4, use3d, solvent);
		        	                        					}
		        	                                        	predictioncount++;
		        	                                        	spheres+=resulth[4];
		        		                                    	if(!donetocsy.contains(resultc[1]+";"+resulth[1])) {
		        		                                    		//System.out.format(Locale.US, "%3d   %3d%8.2f%8.2f\n", i+1, mol.getAtomNumber(away3atom)+1, resultc[1], resulth[1]);
		        		                                    		//System.out.println(resultc[1]+" "+resulth[1]);
		        		                                    		fos.write(new String(resultc[1]+","+resulth[1]+",t\n").getBytes());
		        		                                    		if(debug)
		        		                                    			foshsqctocsy.write(new String(resultc[1]+","+resulth[1]+"\n").getBytes());
		        		                                        	donetocsy.add(resultc[1]+";"+resulth[1]);
		        		                                    	}
		        	                        				}
		        	                        			}
		        	                        			break;
		        	                        		}
		                                		}
		                                	}
		                                	break;
		                        		}
	                        		}
	                        	}
                        	}
                        }
                	}
                    if(debug) {
                    	foshsqc.write(new String("Quality "+(spheres/predictioncount)+"\n").getBytes());
                    	if(doHmbc)
                    		foshmbc.write(new String("Quality "+(spheres/predictioncount)+"\n").getBytes());
                    	if(doHsqcTocsy)
                    		foshsqctocsy.write(new String("Quality "+(spheres/predictioncount)+"\n").getBytes());
                    }
                	fos.write(new String("/\n").getBytes());
                	k++;
                }
                fos.close();
                if(debug) {
                	foshsqc.close();
                	if(doHmbc)
                		foshmbc.close();
                	if(doHsqcTocsy)
                		foshsqctocsy.close();
                }
                smilesreader.close();
            	java.lang.System.exit(0);
        }
        
        private static float[] shift(JSONArray predsmol, int atomid) {
        	for(int i=0;i<predsmol.size();i++) {
        		if(((Long)((JSONObject)predsmol.get(i)).get("atom_idx")).intValue()==atomid)
        			return new float[] {0,((Double)((JSONObject)predsmol.get(i)).get("pred_mu")).floatValue(),0,0,0};
        	}
			return null;
		}
        
		public static void addAndPlaceHydrogens(IAtomContainer container) throws IOException, ClassNotFoundException, CDKException {
        	addExplicitHydrogensToSatisfyValency(container);
            for (int i = 0; i < container.getAtomCount(); i++) {
              IAtom a = container.getAtom(i);
              if (container.getConnectedAtomsList(a).size() == 4) {
                int up = 0;
                int down = 0;
                int hs = 0;
                IAtom h = null;
                for (int k = 0; k < 4; k++) {
                  if (container.getBond(a, (IAtom)container.getConnectedAtomsList(a).get(k)).getStereo() == IBond.Stereo.UP) {
                    up++;
                  }
                  if (container.getBond(a, (IAtom)container.getConnectedAtomsList(a).get(k)).getStereo() == IBond.Stereo.DOWN) {
                    down++;
                  }
                  if (container.getBond(a, (IAtom)container.getConnectedAtomsList(a).get(k)).getStereo() == IBond.Stereo.NONE && ((IAtom)container.getConnectedAtomsList(a).get(k)).getSymbol().equals("H")) {
                    h = (IAtom)container.getConnectedAtomsList(a).get(k);
                    hs++;
                  }
                }
                if (up == 0 && down == 1 && h != null && hs == 1) {
                  container.getBond(a, h).setStereo(IBond.Stereo.UP);
                }
                if (up == 1 && down == 0 && h != null && hs == 1) {
                  container.getBond(a, h).setStereo(IBond.Stereo.DOWN);
                }
              }
            }
            double bondLength = GeometryTools.getBondLengthAverage(container);
            new HydrogenPlacer().placeHydrogens2D(container, bondLength);
        }
        
        public static void addExplicitHydrogensToSatisfyValency(IAtomContainer mol) throws CDKException{
            addImplicitHydrogensToSatisfyValency(mol);
        	AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
        }
        
        public static void addImplicitHydrogensToSatisfyValency(IAtomContainer mol) throws CDKException{
        	CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(mol.getBuilder());    	
            CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(mol.getBuilder());
        	for (IAtom atom : mol.atoms()) {
        		IAtomType type = matcher.findMatchingAtomType(mol, atom);
        		if(type!=null){
        		    AtomTypeManipulator.configure(atom, type);
        		    hAdder.addImplicitHydrogens(mol, atom);
        		}
        	}
          }
}
