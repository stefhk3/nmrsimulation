package uk.ac.dmu.simulate;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.Properties;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.iterator.IteratingSMILESReader;
import org.openscience.cdk.layout.StructureDiagramGenerator;

/**
 * Converts a smiles file to an sdf. Hydrogens are added. Explicit hydrogens are removed first to have all hydrogens at the end of the atom list.
 *
 */
public class Convert {

	public static void main(String[] args) throws CDKException, FileNotFoundException, IOException, ClassNotFoundException {
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
		if(props.get("usedeeplearning").equals("true")) {
	        IteratingSMILESReader smilesreader = new IteratingSMILESReader(new FileInputStream(new File(projectdir+props.getProperty("msmsinput"))), DefaultChemObjectBuilder.getInstance());
	        File sdf=new File(projectdir+props.getProperty("msmsinput").substring(0,props.getProperty("msmsinput").length()-4)+".sdf");
	        FileOutputStream fos=new FileOutputStream(sdf);
	        SDFWriter sdfwriter=new SDFWriter(fos);
	        while(smilesreader.hasNext()) {
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
	        	Simulate.addAndPlaceHydrogens(mol);
	        	sdfwriter.write(mol);
	        }
	        sdfwriter.close();
	        smilesreader.close();
		}
		System.out.print(props.get("usedeeplearning"));
	}

}
