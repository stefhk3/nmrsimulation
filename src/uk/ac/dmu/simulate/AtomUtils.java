/*  This file is part of nmrshiftdb2, an NMR database package.
 *  Copyright (C) 2010-2015 Stefan Kuhn
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package uk.ac.dmu.simulate;

import java.io.IOException;
import java.io.StringReader;
import java.math.BigInteger;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.config.Isotopes;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.geometry.BondTools;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.graph.invariant.MorganNumbersTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.io.SMILESReader;
import org.openscience.cdk.layout.HydrogenPlacer;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.qsar.descriptors.atomic.BondsToAtomDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.EffectiveAtomPolarizabilityDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.IsProtonInAromaticSystemDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.IsProtonInConjugatedPiSystemDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.PartialPiChargeDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.PartialSigmaChargeDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.PeriodicTablePositionDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.PiElectronegativityDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.RDFProtonDescriptor_G3R;
import org.openscience.cdk.qsar.descriptors.atomic.RDFProtonDescriptor_GDR;
import org.openscience.cdk.qsar.descriptors.atomic.RDFProtonDescriptor_GHR;
import org.openscience.cdk.qsar.descriptors.atomic.RDFProtonDescriptor_GHR_topol;
import org.openscience.cdk.qsar.descriptors.atomic.RDFProtonDescriptor_GSR;
import org.openscience.cdk.qsar.descriptors.atomic.SigmaElectronegativityDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.VdWRadiusDescriptor;
import org.openscience.cdk.qsar.descriptors.atompair.PiContactDetectionDescriptor;
import org.openscience.cdk.smiles.FixBondOrdersTool;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.HOSECodeGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import org.w3c.dom.Document;

/**
 *Contains some utils for Molecules/Atoms
 *
 * @author     shk3
 * @created    12. April 2002
 */
public class AtomUtils {

  static String rowstext = "";
  static int rownumber = 0;
  static int colsnumber = 0;
  static VdWRadiusDescriptor vdwrdes = null;
  static SigmaElectronegativityDescriptor sedes = null;
  static PiElectronegativityDescriptor pedes = null;
  static BondsToAtomDescriptor btades = new BondsToAtomDescriptor();
  static PeriodicTablePositionDescriptor ptpdes = new PeriodicTablePositionDescriptor();
  final static Boolean FALSE = new Boolean(false);
  static HOSECodeGenerator hcg = null;
  static RDFProtonDescriptor_G3R rdfg3r = null;
  static RDFProtonDescriptor_GDR rdfgdr = null;
  static RDFProtonDescriptor_GHR rdfghr = null;
  static RDFProtonDescriptor_GSR rdfgsr = null;
  static RDFProtonDescriptor_GHR_topol rdfghrtopol = null;
  static IsProtonInConjugatedPiSystemDescriptor ipicpsdes = null;
  static IsProtonInAromaticSystemDescriptor ipiasdes = null;
  static EffectiveAtomPolarizabilityDescriptor epdes = null;
  static PartialPiChargeDescriptor ppichargedes=null;
  static PartialSigmaChargeDescriptor psigmachargedes=null;
  static PiContactDetectionDescriptor picondes = new PiContactDetectionDescriptor();
  private static Document normalizerdoc = null;
  private static IAtomContainer oldmol = null;
  private static IRingSet allrings = null;
  private static Integer[] maxparameters = {new Integer(20)};



  /**
   *Constructor for the AtomUtils object (private since only static methods)
   */
  private AtomUtils() { }


  /**
   *Gets a BigIntger value representation of the Fingerprint bitset
   *
   * @param  bs   Fingerprint bitset
   * @param  num  A number between 63 and 1023 representing one of the 16 sectors of the Fingerprint
   * @return      The BigInteger representations of the specified Fingerprint sector
   */
  public static BigInteger getBigIntegerValue(BitSet bs, int num) {
    BigInteger bi = new BigInteger("0");
    for (int i = 0; i < 64; i++) {
      if (bs.get(i + (num * 64))) {
    	  bi = bi.add(new BigInteger("2").pow(i));
      }
    }
    return bi;
  }


  /**
   *  Gets the proton Class of a proton according to aires-do-sousa/gasteiger system.
   *
   * @param  cdkmol         The molecule to analyze
   * @param  atomnumber     The atom to classify
   * @param  arf            An existing instace of allringsfinder, to save time.
   * @return                The proton class (1-4)
   * @exception  Exception  Description of Exception
   */
  public static int getProtonClass(IAtomContainer cdkmol, int atomnumber, IRingSet ringSetL) {
    if (cdkmol != oldmol) {
      allrings=ringSetL;
    }
    int protonclass = 4;
    boolean class2 = false;
    for (int l = 0; l < cdkmol.getBondCount(); l++) {
      if (cdkmol.getBond(l).getAtom(0) == cdkmol.getConnectedAtomsList(cdkmol.getAtom(atomnumber)).get(0) || cdkmol.getBond(l).getAtom(1) == cdkmol.getConnectedAtomsList(cdkmol.getAtom(atomnumber)).get(0)) {
        if (cdkmol.getBond(l).getOrder() == IBond.Order.DOUBLE) {
          class2 = true;
        }
      }
    }
    if (((IAtom)cdkmol.getConnectedAtomsList(cdkmol.getAtom(atomnumber)).get(0)).getFlag(CDKConstants.ISAROMATIC)) {
      protonclass = 1;
    } else if (allrings.getRings((IAtom)cdkmol.getConnectedAtomsList(cdkmol.getAtom(atomnumber)).get(0)).getAtomContainerCount() > 0) {
      protonclass = 3;
    } else if (class2) {
      protonclass = 2;
    }
    return protonclass;
  }

  /**
   *  Gets the carbon Class of a proton according to aires-do-sousa/gasteiger system.
   *
   * @param  cdkmol         The molecule to analyze
   * @param  atomnumber     The atom to classify
   * @param  arf            An existing instace of allringsfinder, to save time.
   * @return                The proton class (1-4)
   * @exception  Exception  Description of Exception
   */
  public static int getCarbonClass(IAtomContainer cdkmol, int atomnumber, IRingSet ringSetL) throws Exception {
    if (cdkmol != oldmol) {
      allrings=ringSetL;
    }
    int protonclass = 4;
    boolean class2 = false;
    for (int l = 0; l < cdkmol.getBondCount(); l++) {
      if (cdkmol.getBond(l).getAtom(0) == cdkmol.getAtom(atomnumber) || cdkmol.getBond(l).getAtom(1) == cdkmol.getAtom(atomnumber)) {
        if (cdkmol.getBond(l).getOrder() == IBond.Order.DOUBLE) {
          class2 = true;
        }
      }
    }
    if (cdkmol.getAtom(atomnumber).getFlag(CDKConstants.ISAROMATIC)) {
      protonclass = 1;
    } else if (allrings.getRings(cdkmol.getAtom(atomnumber)).getAtomContainerCount() > 0) {
      protonclass = 3;
    } else if (class2) {
      protonclass = 2;
    }
    return protonclass;
  }

  /**
   *  Gets the allRings attribute, set with last run of getProtonClass.
   *
   * @return    The allRings value
   */
  public static IRingSet getAllRings() {
    return allrings;
  }


  private static void initminusone(float[] floats){
	  for(int i=0;i<floats.length;i++){
		  floats[i]=-1;
	  }
  }

  /**
   * Parses a fuzzy formula search criteria string.
   *
   * @param  criteria                    The criteria to analyze.
   * @param  checkElements               Tells if element name should be checked for validity.
   * @param  checkmultis                 Description of Parameter
   * @return                             A map with elements names as keys and an array of two strings as value which are lower and upper limit.
   * @exception  NmrshiftdbException     Criteria not correct (see message).
   * @exception  IOException             From IsotopeFactory.
   * @exception  ClassNotFoundException  From IsotopeFactory.
   */
  public static Map<String, String[]> parseFuzzyCriteria(String criteria, boolean checkElements, boolean checkmultis) throws NmrshiftdbException, IOException, ClassNotFoundException {
    if (Character.isLetter(criteria.charAt(criteria.length() - 1))) {
      criteria = criteria + "1";
    }
    if (!Character.isLetter(criteria.charAt(0))) {
      throw new NmrshiftdbException("Syntax error");
    }
    boolean addingChars = true;
    HashMap<String,String[]> result = new HashMap<String,String[]>();
    StringBuffer chars = new StringBuffer();
    StringBuffer ranges = new StringBuffer();
    for (int i = 0; i < criteria.length(); i++) {
      if (addingChars && i + 1 < criteria.length() && (!Character.isLetter(criteria.charAt(i + 1)) || Character.isUpperCase(criteria.charAt(i + 1)))) {
        addingChars = false;
        chars.append(criteria.charAt(i));
      }
      if (addingChars) {
        chars.append(criteria.charAt(i));
        if (i + 1 == criteria.length()) {
          addingChars = false;
        }
      } else {
        if (!Character.isLetter(criteria.charAt(i))) {
          ranges.append(criteria.charAt(i));
        }
        if (i == criteria.length() - 1 || (i + 1 < criteria.length() && Character.isLetter(criteria.charAt(i + 1)))) {
          addingChars = true;
          if (!chars.toString().equals("")) {
            if (ranges.toString().equals("")) {
              ranges.append("1");
            }
            if (checkElements && Isotopes.getInstance().getElement(chars.toString()) == null) {
              throw new NmrshiftdbException("Invalid element name");
            }
            String lowerLimit = "";
            String upperLimit = "";
            if (ranges.toString().indexOf("-") == -1) {
              if (ranges.toString().equals("*")) {
                lowerLimit = "0";
                upperLimit = "" + Integer.MAX_VALUE;
              } else {
                try {
                  lowerLimit = "" + Integer.parseInt(ranges.toString());
                  upperLimit = lowerLimit;
                } catch (NumberFormatException ex) {
                  throw new NmrshiftdbException("One of the ranges contains no -, but is not a figure as well");
                }
              }
            } else {
              try {
                lowerLimit = "" + Integer.parseInt(ranges.toString().substring(0, ranges.toString().indexOf("-")));
                upperLimit = "" + Integer.parseInt(ranges.toString().substring(ranges.toString().indexOf("-") + 1, ranges.toString().length()));
              } catch (NumberFormatException ex) {
                throw new NmrshiftdbException("One of the ranges contains a -, but before or after it is not a figure!");
              }
            }
            if (result.get(chars.toString()) != null) {
            	result.put(chars.toString(), new String[]{new String(""+(Integer.parseInt(lowerLimit)+Integer.parseInt(((String[])result.get(chars.toString()))[0]))), new String(""+(upperLimit+Integer.parseInt(((String[])result.get(chars.toString()))[1])))});
            }else{
            	result.put(chars.toString(), new String[]{lowerLimit, upperLimit});
            }
          } else {
            throw new NmrshiftdbException("Syntax error");
          }
          chars = new StringBuffer();
          ranges = new StringBuffer();
        }
      }
    }
    if (!checkElements && checkmultis) {
      if (result.get("S") == null && result.get("s") == null) {
        result.put("S", new String[]{"0", "0"});
      }
      if (result.get("D") == null && result.get("d") == null) {
        result.put("D", new String[]{"0", "0"});
      }
      if (result.get("T") == null && result.get("t") == null) {
        result.put("T", new String[]{"0", "0"});
      }
      if (result.get("Q") == null && result.get("q") == null) {
        result.put("Q", new String[]{"0", "0"});
      }
    }
    return (result);
  }
  
  /**
   * Adds implicit hydrogens to a molecule to satisfy valency rules.
   * 
   * @param mol The molecule to saturate.
   * @throws CDKException
   */
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


  /**
   * Adds explicit hydrogens to a molecule to satisfy valency rules (no layout done).
   * 
   * @param mol The molecule to saturate.
   * @throws CDKException
   */
  public static void addExplicitHydrogensToSatisfyValency(IAtomContainer mol) throws CDKException{
    addImplicitHydrogensToSatisfyValency(mol);
	AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
  }

  /**
   *  Adds and places hydrogens to a molecule. Also makes up/down bonds if there is 
   *  only one up/down bond and one h an a stereocenter.
   *
   * @param  container                   The molecule to saturate.
   * @exception  IOException             Description of Exception.
   * @exception  ClassNotFoundException  Description of Exception.
   * @exception  CDKException            Description of Exception
   */
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


  /**
   *  Sorts atoms like Cs-other heavy atoms-Hs.
   *
   * @param  mol            The molecule to rearrange.
   * @return                The new molecule as mdl file.
   */
  public static IAtomContainer rearangeAtoms(IAtomContainer mol) {
    Iterator<IAtom> atomsold = mol.atoms().iterator();
    IAtom[] atomsnew = new IAtom[mol.getAtomCount()];
    int k = 0;
    while(atomsold.hasNext()) {
      IAtom atom=atomsold.next();
      if (atom.getSymbol().equals("C")) {
        atomsnew[k++] = atom;
      }
    }
    atomsold = mol.atoms().iterator();
    while(atomsold.hasNext()) {
      IAtom atom=(IAtom)atomsold.next();
      if (!atom.getSymbol().equals("C") && !atom.getSymbol().equals("H")) {
        atomsnew[k++] = atom;
      }
    }
    atomsold = mol.atoms().iterator();
    while(atomsold.hasNext()) {
      IAtom atom=(IAtom)atomsold.next();
      if (atom.getSymbol().equals("H")) {
        atomsnew[k++] = atom;
      }
    }
    mol.setAtoms(atomsnew);
    return (mol);
  }

  /**
   * Converts in int into an IBond.Order value. 1 equals single, 2 equals double and so on.
   * 
   * @param order The value to translate.
   * @return The order value.
   */
  public static IBond.Order getOrderFromInt(int order){
	  return order==1 ? IBond.Order.SINGLE : order==2 ? IBond.Order.DOUBLE : order ==3 ? IBond.Order.TRIPLE : IBond.Order.QUADRUPLE;
	  
  }


  /**
   * Translates an int into an IBond.Stereo value. It is the reverse of IBond.Stero.ordinal(), so
   * getStereo(x).ordinal()==x. 
   * 
   * @param stereo The int value to translate
   * @return The stereo value.
   * @throws NmrshiftdbException There is no stereo representation for that int. 
   */
  public static IBond.Stereo getStereoFromInt(int stereo) throws NmrshiftdbException {
        if(stereo==IBond.Stereo.NONE.ordinal())
                return IBond.Stereo.NONE;
        else if(stereo==IBond.Stereo.UP.ordinal())
                return IBond.Stereo.UP;
        else if(stereo==IBond.Stereo.UP_INVERTED.ordinal())
                return IBond.Stereo.UP_INVERTED;
        else if(stereo==IBond.Stereo.DOWN.ordinal())
                return IBond.Stereo.DOWN;
        else if(stereo==IBond.Stereo.DOWN_INVERTED.ordinal())
                return IBond.Stereo.DOWN_INVERTED;
        else if(stereo==IBond.Stereo.UP_OR_DOWN.ordinal())
                return IBond.Stereo.UP_OR_DOWN;
        else if(stereo==IBond.Stereo.UP_OR_DOWN_INVERTED.ordinal())
            return IBond.Stereo.UP_OR_DOWN_INVERTED;
        else if(stereo==IBond.Stereo.E_OR_Z.ordinal())
            return IBond.Stereo.E_OR_Z;
        else if(stereo==IBond.Stereo.E.ordinal())
            return IBond.Stereo.E;
        else if(stereo==IBond.Stereo.Z.ordinal())
            return IBond.Stereo.Z;
        else if(stereo==IBond.Stereo.E_Z_BY_COORDINATES.ordinal())
            return IBond.Stereo.E_Z_BY_COORDINATES;
        else throw new NmrshiftdbException("No such stereo value: "+stereo); 
  }


  	/**
	 * Gets an atom from a moleucle by ID (not atom number), null if none found. 
	 * 
	 * @param molWithH The molecule.
	 * @param string   The ID.
	 * @return The atom found or null.
	 */
	public static IAtom getAtomById(IAtomContainer molWithH, String string) {
	 for(int i=0;i<molWithH.getAtomCount();i++){
	  if(molWithH.getAtom(i).getID().equals(string))
		  return molWithH.getAtom(i);
	 }
	 return null;
	}
	
	/**
	 * Reads a molecule from smles and lays it out using cdk StructureDiagramGenerator.
	 * 
	 * @param smiles The smiles to use.
	 * @return The molecule.
	 * @throws CDKException
	 * @throws IOException 
	 */
	public static IAtomContainer readFromSmilesAndLayout(String smiles) throws CDKException, IOException{
		SMILESReader reader = new SMILESReader(new StringReader(smiles));
        AtomContainerSet som = (AtomContainerSet)reader.read(new AtomContainerSet());
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(som.getAtomContainer(0));
        FixBondOrdersTool dbst = new FixBondOrdersTool();
        IAtomContainer molecule = dbst.kekuliseAromaticRings(som.getAtomContainer(0));
        StructureDiagramGenerator sdg = new StructureDiagramGenerator();
        sdg.setMolecule(molecule);
        sdg.generateCoordinates();
        reader.close();
        return sdg.getMolecule();
	}

	
	//Improved version of cdk BondTools.makeUpDownBonds
    public static void  makeUpDownBonds(IAtomContainer container){
	    for (int i = 0; i < container.getAtomCount(); i++) {
	        IAtom a = container.getAtom(i);
	        if (container.getConnectedAtomsList(a).size() == 4) {
	          int up = 0;
	          int down = 0;
	          int hs = 0;
	          IAtom h = null;
	          for (int k = 0; k < 4; k++) {
	        	  IAtom conAtom = container.getConnectedAtomsList(a).get(k);
	        	  IBond conBond = container.getBond(a,conAtom);
	        	  IBond.Stereo stereo = conBond.getStereo();
	        	  if (stereo == IBond.Stereo.UP && conBond.getAtom(1)==conAtom || stereo == IBond.Stereo.UP_INVERTED && conBond.getAtom(0)==conAtom) {
	        		  up++;
	        	  }
	        	  else if (stereo == IBond.Stereo.DOWN && conBond.getAtom(1)==conAtom || stereo == IBond.Stereo.DOWN_INVERTED && conBond.getAtom(0)==conAtom) {
	        		  down++;
	        	  }
	        	  else if (stereo == IBond.Stereo.NONE && conAtom.getSymbol().equals("H")) {
	        		  h = conAtom;
	        		  hs++;
	        	  } else {
	        		  h = null;
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
	}

    
    public static int isTetrahedral(IAtomContainer container, IAtom atom, boolean strict)
	{
		List<IAtom> atoms = container.getConnectedAtomsList(atom);
		if (atoms.size() != 4)
		{
			return (0);
		}
		java.util.List<IBond> bonds = container.getConnectedBondsList(atom);
        int up = 0;
		int down = 0;
        for (IBond bond : bonds) {
            if (bond.getStereo() == IBond.Stereo.NONE ||
            	bond.getStereo() == CDKConstants.UNSET) {
            } else
            if ((bond.getStereo() == IBond.Stereo.UP && bond.getAtom(0)==atom) || (bond.getStereo() == IBond.Stereo.UP_INVERTED && bond.getAtom(1)==atom)) {
                up++;
            } else
            if ((bond.getStereo() == IBond.Stereo.DOWN && bond.getAtom(0)==atom) || (bond.getStereo() == IBond.Stereo.DOWN_INVERTED && bond.getAtom(1)==atom)) {
                down++;
            }
        }
		if (up == 1 && down == 1)
		{
			return 1;
		}
		if (up == 2 && down == 2)
		{
			if (BondTools.stereosAreOpposite(container, atom))
			{
				return 2;
			}
			return 0;
		}
		if (up == 1 && down == 0 && !strict)
		{
			return 3;
		}
		if (down == 1 && up == 0 && !strict)
		{
			return 4;
		}
		if (down == 2 && up == 1 && !strict)
		{
			return 5;
		}
		if (down == 1 && up == 2 && !strict)
		{
			return 6;
		}
		return 0;
	}
    
    
    /**
     *  Tells if a certain bond is center of a valid double bond configuration.
     *  Copy of cdk BondTools method, with improvements
     *
     * @param  container  The atomcontainer.
     * @param  bond       The bond.
     * @param  sssr       The smallestSetOfSmallestRings of the container.
     * @return            true=is a potential configuration, false=is not.
     */
    public static boolean isValidDoubleBondConfiguration(IAtomContainer container, IBond bond, IRingSet sssr) {
    	if(bond.getFlag(CDKConstants.ISAROMATIC) || (sssr!=null && containsInRingSmallerThanEight(sssr, bond)))
    		return false;
      List<IAtom> connectedAtoms = container.getConnectedAtomsList(bond.getAtom(0));
      IAtom from = null;
        for (IAtom connectedAtom : connectedAtoms) {
            if (connectedAtom != bond.getAtom(1)) {
                from = connectedAtom;
            }
        }
      boolean[] array = new boolean[container.getBondCount()];
      for (int i = 0; i < array.length; i++) {
        array[i] = true;
      }
      if (isStartOfDoubleBond(container, bond.getAtom(0), from, array) && isEndOfDoubleBond(container, bond.getAtom(1), bond.getAtom(0), array) && !bond.getFlag(CDKConstants.ISAROMATIC)) {
        return (true);
      } else {
        return (false);
      }
    }

    
    /**
     * Tells if this ring set contains this bond in a ring of a size smaller eight atoms.
     * 
     * @param sssr The ring set.
     * @param bond The bond.
     * @return true=is in a ring smaller eight, false=is not.
     */
    private static boolean containsInRingSmallerThanEight(IRingSet sssr,
			IBond bond) {
		for(IAtomContainer ring : sssr.atomContainers()){
			if(ring.getAtomCount()<8 && ring.contains(bond))
				return true;
		}
		return false;
	}


	/**
     *  Says if an atom is the end of a double bond configuration
     *  Copy of cdk BondTools method, with improvements
     *
     * @param  atom                     The atom which is the end of configuration
     * @param  container                The atomContainer the atom is in
     * @param  parent                   The atom we came from
     * @param  doubleBondConfiguration  The array indicating where double bond
     *      configurations are specified (this method ensures that there is
     *      actually the possibility of a double bond configuration)
     * @return                          false=is not end of configuration, true=is
     */
    private static boolean isEndOfDoubleBond(IAtomContainer container, IAtom atom, IAtom parent, boolean[] doubleBondConfiguration) {
      if (container.getBondNumber(atom, parent) == -1 || doubleBondConfiguration.length <= container.getBondNumber(atom, parent) || !doubleBondConfiguration[container.getBondNumber(atom, parent)]) {
        return false;
      }

        int hcount;
        if (atom.getImplicitHydrogenCount() == CDKConstants.UNSET) hcount = 0;
        else hcount = atom.getImplicitHydrogenCount();

      int lengthAtom = container.getConnectedAtomsList(atom).size() + hcount;

         if (parent.getImplicitHydrogenCount() == CDKConstants.UNSET) hcount = 0;
        else hcount = parent.getImplicitHydrogenCount();

      int lengthParent = container.getConnectedAtomsList(parent).size() + hcount;
        
      if (container.getBond(atom, parent) != null) {
        if (container.getBond(atom, parent).getOrder() == Order.DOUBLE && (lengthAtom == 3 || (lengthAtom == 2 && atom.getSymbol().equals("N"))) && (lengthParent == 3 || (lengthParent == 2 && parent.getSymbol().equals("N")))) {
          List<IAtom> atoms = container.getConnectedAtomsList(atom);
          IAtom one = null;
          IAtom two = null;
            for (IAtom conAtom : atoms) {
                if (conAtom != parent && one == null) {
                    one = conAtom;
                } else if (conAtom != parent && one != null) {
                    two = conAtom;
                }
            }
          String[] morgannumbers = MorganNumbersTools.getMorganNumbersWithElementSymbol(container);
          if ((one != null && two == null && atom.getSymbol().equals("N") && Math.abs(BondTools.giveAngleBothMethods(parent, atom, one, true)) > Math.PI / 10) || (!atom.getSymbol().equals("N") && one != null && two != null && !morgannumbers[container.getAtomNumber(one)].equals(morgannumbers[container.getAtomNumber(two)]))) {
            return (true);
          } else {
            return (false);
          }
        }
      }
      return (false);
    }


    /**
     *  Says if an atom is the start of a double bond configuration
     *  Copy of cdk BondTools method, with improvements
     *
     * @param  a                        The atom which is the start of configuration
     * @param  container                The atomContainer the atom is in
     * @param  parent                   The atom we came from
     * @param  doubleBondConfiguration  The array indicating where double bond
     *      configurations are specified (this method ensures that there is
     *      actually the possibility of a double bond configuration)
     * @return                          false=is not start of configuration, true=is
     */
    private static boolean isStartOfDoubleBond(IAtomContainer container, IAtom a, IAtom parent, boolean[] doubleBondConfiguration) {
        int hcount;
        if (a.getImplicitHydrogenCount() == CDKConstants.UNSET) hcount = 0;
        else hcount = a.getImplicitHydrogenCount();

      int lengthAtom = container.getConnectedAtomsList(a).size() + hcount;
        
      if (lengthAtom != 3 && (lengthAtom != 2 && !(a.getSymbol().equals("N")))) {
        return (false);
      }
      List<IAtom> atoms = container.getConnectedAtomsList(a);
      IAtom one = null;
      IAtom two = null;
      boolean doubleBond = false;
      IAtom nextAtom = null;
        for (IAtom atom : atoms) {
            if (atom != parent && container.getBond(atom, a).getOrder() == Order.DOUBLE && isEndOfDoubleBond(container, atom, a, doubleBondConfiguration)) {
                doubleBond = true;
                nextAtom = atom;
            }
            if (atom != nextAtom && one == null) {
                one = atom;
            } else if (atom != nextAtom && one != null) {
                two = atom;
            }
        }
      String[] morgannumbers = MorganNumbersTools.getMorganNumbersWithElementSymbol(container);
      if (one != null && ((!a.getSymbol().equals("N") && two != null && !morgannumbers[container.getAtomNumber(one)].equals(morgannumbers[container.getAtomNumber(two)]) && doubleBond && doubleBondConfiguration[container.getBondNumber(a, nextAtom)]) || (doubleBond && a.getSymbol().equals("N") && Math.abs(BondTools.giveAngleBothMethods(nextAtom, a, parent, true)) > Math.PI / 10))) {
      	return (true);
      } else {
      	return (false);
      }
    }


	public static int getHcount(IAtomContainer mol, IAtom atom) {
		List<IAtom> atoms = mol.getConnectedAtomsList(atom);
		int NumberOfHs=0;
        for (int k = 0; k < atoms.size(); k++) {
          if (atoms.get(k).getSymbol().equals("H")) {
            NumberOfHs++;
          }
        }		
		return NumberOfHs;
	}
}

