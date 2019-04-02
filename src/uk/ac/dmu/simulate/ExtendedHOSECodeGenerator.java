/*  This file is part of nmrshiftdb2, an NMR database package.
 *  Copyright (C) 1997-2007  Christoph Steinbeck <steinbeck@users.sf.net>
 *                     2008  Egon Willighagen <egonw@users.sf.net>
 *  2013 Stefan Kuhn <shk3@users.sf.net>
 *
 *  This code was written by Stefan Kuhn using the class org.openscience.
 *  cdk.tools.HOSECodeGenerator from the Chemistry Development Kit 
 *  (sf.net/projects/cdk) rev. [2fc6b6] (if I got git right). This new 
 *  class is licensed with all of nmrshiftdb2 under AGPL v. 3, the original 
 *  software is of course still available under lGPL.
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
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.TreeMap;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.config.IsotopeFactory;
import org.openscience.cdk.config.Isotopes;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.geometry.BondTools;
import org.openscience.cdk.graph.invariant.CanonicalLabeler;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IIsotope;
import org.openscience.cdk.interfaces.IRing;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.SSSRFinder;
import org.openscience.cdk.smiles.InvPair;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 * Generates HOSE codes {@cdk.cite BRE78}.
 * IMPORTANT: Your molecule must contain implicit or explicit hydrogens
 * for this method to work properly.
 *
 * @author     steinbeck
 * @cdk.githash
 * @cdk.keyword    HOSE code, spherical atom search
 * @cdk.created    2002-05-10
 * @cdk.module     standard
 */
public class ExtendedHOSECodeGenerator implements java.io.Serializable
{

	private static ILoggingTool logger =
	    LoggingToolFactory.createLoggingTool(ExtendedHOSECodeGenerator.class);
	
    private static final long serialVersionUID = -4353471818831864513L;
    
    /**
	 *  Container for the nodes in a sphere.
	 */
	protected List<TreeNode> sphereNodes = null;
    protected List<IAtom> sphereNodesWithAtoms = null;
    
	/**
	 *  Container for the node in the next sphere Assembled in a recursive method
	 *  and then passed to the next recursion to become "sphereNodes".
	 */
	protected List<TreeNode> nextSphereNodes = null;
	/**
	 *  Counter for the sphere in which we currently work.
	 */
	protected int sphere = 0;
	/**
	 *  How many spheres are we supposed inspect.
	 */
	protected int maxSphere = 0;

	/**
	 *  Here we store the spheres that we assemble, in order to parse them into a
	 *  code later.
	 */
	protected List<TreeNode>[] spheres = null;
	protected List<IAtom>[] spheresWithAtoms = null;

	/**
	 *  The HOSECode string that we assemble
	 */
	protected StringBuffer HOSECode = null;

	/**
	 *  The molecular structure on which we work
	 */
	protected IAtomContainer atomContainer;

	/**
	 *  Delimiters used to separate spheres in the output string. Bremser uses the
	 *  sequence"(//)" for the first four spheres.
	 */
	protected String[] sphereDelimiters =
			{
			"(", "/", "/", ")", "/", "/", "/", "/", "/", "/", "/", "/"
			};
	/**
	 *  The bond symbols used for bond orders "single", "double", "triple" and
	 *  "aromatic"
	 */
	protected String bondSymbols[] =
			{
			"", "", "=", "%", "*"
			};
			
	protected String centerCode = null;			
			
	public TreeNode rootNode = null;
			
			
	boolean debug = false;
	
	private IAtomContainer acold=null;
	private IRingSet soar=null;

	/**
	 *  The rank order for the given element symbols.
	 */

	static final String[] rankedSymbols =
			{
			"C", "O", "N", "S", "P", "Si", "B", "F", "Cl", "Br", ";", "I", "#", "&", ","
			};

	/**
	 *  The ranking values to be used for the symbols above.
	 */
	static final int[] symbolRankings =
			{
			9000, 8900, 8800, 8700,
			8600, 8500, 8400, 8300, 8200, 8100, 8000, 7900, 1200, 1100, 1000
			};

	/**
	 *  The bond rankings to be used for the four bond order possibilities.
	 */

	static final int[] bondRankings =
			{
			0, 0, 200000, 300000, 100000
			};


	/**
	 *  Constructor for the HOSECodeGenerator.
	 */
	public ExtendedHOSECodeGenerator()
	{
		sphereNodes = new ArrayList<TreeNode>();
		sphereNodesWithAtoms = new ArrayList<IAtom>();
		nextSphereNodes = new ArrayList<TreeNode>();
		HOSECode = new StringBuffer();
	}
  
	private IsotopeFactory isotopeFac = null;
	
	private void ensureIsotopeFactory() throws CDKException {
		if (isotopeFac == null) {
			try {
	            isotopeFac = Isotopes.getInstance();
            } catch (IOException e) {
	            throw new CDKException("Could not instantiate the IsotopeFactory: " + e.getMessage(), e);
            }
		}
	}
  
	/**
	 *  This method is intended to be used to get the atoms around an atom in spheres. It is not used in this class, but is provided for other classes to use.
	 *  It also creates the HOSE code in HOSECode as a side-effect.
	 *  
	 *@param  ac  The {@link IMolecule} with the molecular skeleton in which the root atom resides.
	 *@param  root The root atom for which to produce the spheres.
	 *@param  noOfSpheres  The number of spheres to look at.
	 *@param  ringsize  Shall the center code have the ring size in it? Only use if you want to have the hose code later, else say false.
	 *@return An array of {@link List}. The list at i-1 contains the atoms at sphere i as TreeNodes.
	 **/
	public List<IAtom>[] getSpheres(IAtomContainer ac, IAtom root, int noOfSpheres, boolean ringsize) throws CDKException
	{
		ensureIsotopeFactory();
		centerCode = "";
		this.atomContainer = ac;
		maxSphere = noOfSpheres;
		spheres = new List[noOfSpheres + 1];
		spheresWithAtoms = new List[noOfSpheres + 1];
		for (int i = 0; i < ac.getAtomCount(); i++)
		{
			ac.getAtom(i).setFlag(CDKConstants.VISITED, false);
		}
		root.setFlag(CDKConstants.VISITED, true);
		rootNode = new TreeNode(root.getSymbol(), null, root, (double)0, atomContainer.getConnectedBondsCount(root), 0, null, isChiral(root,ac));
		/*
		 *  All we need to observe is how the ranking of substituents
		 *  in the subsequent spheres of the root nodes influences the
		 *  ranking of the first sphere, since the order of a node in a sphere
		 *  depends on the order the preceding node in its branch
		 */
		HOSECode = new StringBuffer();
		createCenterCode(root,ac,ringsize);
		breadthFirstSearch(root, false);
		createCode();
		fillUpSphereDelimiters();
		logger.debug("HOSECodeGenerator -> HOSECode: " + HOSECode.toString());
		return spheresWithAtoms;
	}


	private boolean isChiral(IAtom root, IAtomContainer ac) {
		if(!root.getSymbol().equals("C"))
			return false;
		List<IAtom> conatoms = ac.getConnectedAtomsList(root);
		for(int i=0;i<conatoms.size();i++){
			if(!ac.getBond(root, conatoms.get(i)).getOrder().equals(IBond.Order.SINGLE))
				return false;
		}
		for(int i=0;i<conatoms.size();i++){
			if(ac.getBond(root, conatoms.get(i)).getStereo().equals(IBond.Stereo.DOWN) && ac.getBond(root, conatoms.get(i)).getAtom(0)==root)
				return true;
			if(ac.getBond(root, conatoms.get(i)).getStereo().equals(IBond.Stereo.DOWN_INVERTED) && ac.getBond(root, conatoms.get(i)).getAtom(1)==root)
				return true;
			if(ac.getBond(root, conatoms.get(i)).getStereo().equals(IBond.Stereo.UP) && ac.getBond(root, conatoms.get(i)).getAtom(0)==root)
				return true;
			if(ac.getBond(root, conatoms.get(i)).getStereo().equals(IBond.Stereo.UP_INVERTED) && ac.getBond(root, conatoms.get(i)).getAtom(1)==root)
				return true;
		}
		return false;
	}

	/**
	 * Produces a HOSE code for Atom <code>root</code> in the {@link IAtomContainer} <code>ac</code>. The HOSE
	 * code is produced for the number of spheres given by <code>noOfSpheres</code>.
	 * IMPORTANT: if you want aromaticity to be included in the code, you need
	 * to run the IAtomContainer <code>ac</code> to the {@link CDKHueckelAromaticityDetector} prior to 
	 * using <code>getHOSECode()</code>. This method only gives proper results if the molecule is
	 * fully saturated (if not, the order of the HOSE code might depend on atoms in higher spheres).
	 * This method is known to fail for protons sometimes.
	 * IMPORTANT: Your molecule must contain implicit or explicit hydrogens
	 * for this method to work properly.
	 *
	 * @param  ac  The {@link IAtomContainer} with the molecular skeleton in which the root atom resides
	 * @param  root The root atom for which to produce the HOSE code
	 * @param  noOfSpheres  The number of spheres to look at
	 * @return The HOSECode value
	 * @exception  org.openscience.cdk.exception.CDKException  Thrown if something is wrong
	 */
	public String getHOSECode(IAtomContainer ac, IAtom root, int noOfSpheres) throws CDKException
	{
		return getHOSECode(ac,root,noOfSpheres, false);
	}
	
	
	/**
   * Produces a HOSE code for Atom <code>root</code> in the {@link IAtomContainer} <code>ac</code>. The HOSE
   * code is produced for the number of spheres given by <code>noOfSpheres</code>.
   * IMPORTANT: if you want aromaticity to be included in the code, you need
   * to run the IAtomContainer <code>ac</code> to the {@link CDKHueckelAromaticityDetector} prior to 
   * using <code>getHOSECode()</code>. This method only gives proper results if the molecule is
   * fully saturated (if not, the order of the HOSE code might depend on atoms in higher spheres).
   * This method is known to fail for protons sometimes.
   * IMPORTANT: Your molecule must contain implicit or explicit hydrogens
   * for this method to work properly.
	 *
	 * @param  ac  The IAtomContainer with the molecular skeleton in which the root atom resides
	 * @param  root The root atom for which to produce the HOSE code
	 * @param  noOfSpheres  The number of spheres to look at
	 * @param  ringsize  The size of the ring(s) it is in is included in center atom code
	 * @return The HOSECode value
	 * @exception  org.openscience.cdk.exception.CDKException  Thrown if something is wrong
	 */
	public String getHOSECode(IAtomContainer ac, IAtom root, int noOfSpheres, boolean ringsize) throws CDKException
	{
		ensureIsotopeFactory();
    CanonicalLabeler canLabler = new CanonicalLabeler();
    canLabler.canonLabel(ac);
		centerCode = "";
		this.atomContainer = ac;
		maxSphere = noOfSpheres;
		spheres = new List[noOfSpheres + 1];
		for (int i = 0; i < ac.getAtomCount(); i++)
		{
			ac.getAtom(i).setFlag(CDKConstants.VISITED, false);
		}
		root.setFlag(CDKConstants.VISITED, true);
		rootNode = new TreeNode(root.getSymbol(), null, root, (double)0, atomContainer.getConnectedBondsCount(root), 0, null, isChiral(root,ac));
		/*
		 *  All we need to observe is how the ranking of substituents
		 *  in the subsequent spheres of the root nodes influences the
		 *  ranking of the first sphere, since the order of a node in a sphere
		 *  depends on the order the preceding node in its branch
		 */
		HOSECode = new StringBuffer();
		createCenterCode(root, ac, ringsize);
		breadthFirstSearch(root,true);
		createCode();
		fillUpSphereDelimiters();
		logger.debug("HOSECodeGenerator -> HOSECode: ", HOSECode);
		return HOSECode.toString();
	}

	private void createCenterCode(IAtom root, IAtomContainer ac, boolean ringsize)
	{
		int partnerCount = 0;
		partnerCount = atomContainer.getConnectedBondsCount(root) +
                (root.getImplicitHydrogenCount() == CDKConstants.UNSET ? 0 : root.getImplicitHydrogenCount()); 
		centerCode = root.getSymbol() + "-" + partnerCount + createChargeCode(root)+(ringsize ? getRingcode(root, ac) : "" )+";";
	}
	
	
	private String getRingcode(IAtom root, IAtomContainer ac){
		if(ac!=acold){
			soar=new SSSRFinder(ac).findSSSR();
		}
		boolean[] bool=new boolean[1000];
		StringBuffer sb=new StringBuffer();
		for(int i=0;i<soar.getRings(root).getAtomContainerCount();i++){
			if(((IRing)soar.getRings(root).getAtomContainer(i)).getAtomCount()<bool.length)
				bool[((IRing)soar.getRings(root).getAtomContainer(i)).getAtomCount()]=true;
		}
		for(int i=0;i<bool.length;i++){
			if(bool[i])
				sb.append(i+"");
		}
		if(sb.toString().equals(""))
			return "";
		else
			return "-"+sb.toString();
	}

    private String createChargeCode(IAtom atom) {
        StringBuffer tempCode = new StringBuffer();

        if (atom != null) {
            
            Integer formalCharge = atom.getFormalCharge();
            if (formalCharge == CDKConstants.UNSET) formalCharge = 0;

            if (formalCharge != 0) {

                if (Math.abs(formalCharge) == 1) {
                    if (formalCharge < 0)
                        tempCode.append("-");
                    else
                        tempCode.append("+");
                } else {
                    tempCode.append("'");
                    if (formalCharge > 0)
                        tempCode.append("+");
                    tempCode.append(formalCharge + "'");
                }
            }
        }
        return (tempCode + "");
    }

    /**
	 *  Prepares for a breadth first search within the {@link IAtomContainer}. The actual
	 *  recursion is done in <code>nextSphere()</code>.
	 *
	 *@param  root  The atom at which we start the search
	 *@exception  org.openscience.cdk.exception.CDKException  If something goes wrong.
	 */
	private void breadthFirstSearch(IAtom root,boolean addTreeNode) throws CDKException {
		sphere = 0;
		TreeNode tempNode = null;
		List<IAtom> conAtoms = atomContainer.getConnectedAtomsList(root);
		IAtom atom;
		IBond bond = null;
		sphereNodes.clear();
		sphereNodesWithAtoms.clear();
		for (int i = 0; i < conAtoms.size(); i++){
			try{
				atom = conAtoms.get(i);
				//if(atom.getSymbol().equals("H"))
				//		continue;
				bond = atomContainer.getBond(root, atom);
				/*
				 *  In the first sphere the atoms are labeled with
				 *  their own atom atom as source
				 */
				if (bond.getFlag(CDKConstants.ISAROMATIC))
				{
					tempNode = new TreeNode(atom.getSymbol(), new TreeNode(root.getSymbol(), null, root, (double) 0, 0, (long) 0, null, isChiral(root, atomContainer)), atom, 4, atomContainer.getConnectedBondsCount(atom), 0, bond.getStereo(), isChiral(atom, atomContainer));
				} else
				{
					tempNode = new TreeNode(atom.getSymbol(), new TreeNode(root.getSymbol(), null, root, (double) 0, 0, (long) 0, null, isChiral(root, atomContainer)), atom, bond.getOrder().ordinal()+1, atomContainer.getConnectedBondsCount(atom), 0, bond.getStereo(), isChiral(atom, atomContainer));
				}
				
		        sphereNodes.add(tempNode);
		        if(!addTreeNode)
		        	sphereNodesWithAtoms.add(atom);
		        
//		        rootNode.childs.addElement(tempNode);
				atom.setFlag(CDKConstants.VISITED, true);
			} catch (Exception exc)
			{
				throw new CDKException("Error in HOSECodeGenerator->breadthFirstSearch.", exc);
			}
		}
		Collections.sort(sphereNodes,new TreeNodeComparator());
		nextSphere(sphereNodes);
	}

	/**
	 *  The actual recursion method for our breadth first search. Each node in
	 *  sphereNodes is inspected for its descendants which are then stored in
	 *  <code>nextSphereNodes</code>, which again is passed to the next recursion level of
	 *  <code>nextSphere()</code>.
	 *
	 *@param  sphereNodes The sphereNodes to be inspected
	 *@exception  org.openscience.cdk.exception.CDKException  If something goes wrong
	 */
	private void nextSphere(List<TreeNode> sphereNodes) throws CDKException
	{
		spheres[sphere] = sphereNodes;
		if(spheresWithAtoms!=null)
			spheresWithAtoms[sphere] = sphereNodesWithAtoms;
		/*
		 *  From here we start assembling the next sphere
		 */
        IAtom node = null;
        IAtom toNode = null;
        List<IAtom> conAtoms = null;
		TreeNode treeNode = null;
		nextSphereNodes = new ArrayList<TreeNode>();
		IBond bond = null;
		for (int i = 0; i < sphereNodes.size(); i++)
		{
			treeNode = (TreeNode) sphereNodes.get(i);
			if (!("&;#:,".indexOf(treeNode.symbol) >= 0))
			{
				node = treeNode.atom;
				if(node.getSymbol().equals("H"))
					continue;
				
				conAtoms = atomContainer.getConnectedAtomsList(node);
				if (conAtoms.size() == 1){
					nextSphereNodes.add(new TreeNode(",", treeNode, null, 0, 0, treeNode.score, null, false));
				}else{
					for (int j = 0; j < conAtoms.size(); j++)
					{
						toNode = conAtoms.get(j);
						if (toNode != treeNode.source.atom)
						{
							bond = atomContainer.getBond(node, toNode);
							if (bond.getFlag(CDKConstants.ISAROMATIC))
							{
								nextSphereNodes.add(new TreeNode(toNode.getSymbol(), treeNode, toNode, 4, atomContainer.getConnectedBondsCount(toNode), treeNode.score, bond.getStereo(), isChiral(toNode, atomContainer)));
							} else
							{
								nextSphereNodes.add(new TreeNode(toNode.getSymbol(), treeNode, toNode, bond.getOrder().ordinal()+1, atomContainer.getConnectedBondsCount(toNode), treeNode.score, bond.getStereo(), isChiral(toNode, atomContainer)));
							}
						}
					}
				}
			}
		}
		Collections.sort(nextSphereNodes,new TreeNodeComparator());
		if (sphere < maxSphere)
		{
			sphere++;
			nextSphere(nextSphereNodes);
		}
	}

	public String makeBremserCompliant(String code)
	{
		int sepIndex = code.indexOf(";");
		if (sepIndex >= 0)
		{
			code = code.substring(sepIndex + 1, code.length());	
		}
		return code;
	}
	
	/**
	 *  After recursively having established the spheres and assigning each node an
	 *  appropriate score, we now generate the complete HOSE code.
	 *
	 *@exception  org.openscience.cdk.exception.CDKException  Thrown if something goes wrong
	 */
	private void createCode() throws CDKException {
		List<TreeNode> sphereNodes = null;
		TreeNode tn = null;
		for (int f = 0; f < atomContainer.getAtomCount(); f++)
		{
			atomContainer.getAtom(f).setFlag(CDKConstants.VISITED, false);
		}

		for (int f = 0; f < maxSphere; f++)
		{
			sphereNodes = spheres[maxSphere - f];
			for (int g = 0; g < sphereNodes.size(); g++)
			{
				tn = sphereNodes.get(g);
				if (tn.source != null)
				{
					tn.source.ranking += tn.degree;
				}
				
			}
		}

		for (int f = 0; f < maxSphere; f++)
		{
			sphereNodes = spheres[f];
			calculateNodeScores(sphereNodes);
			sortNodesByScore(sphereNodes);
		}

		for (int f = 0; f < maxSphere; f++)
		{
			sphereNodes = spheres[f];
			for (int g = 0; g < sphereNodes.size() ; g++)
			{
				tn = (TreeNode) sphereNodes.get(g);
				tn.score += tn.ranking;
			}
			sortNodesByScore(sphereNodes);
		}
		for (int f = 0; f < maxSphere; f++)
		{
			sphereNodes = spheres[f];
			List<IAtom> alreadyHandledCenters=new ArrayList<IAtom>();
			for (int g = 0; g < sphereNodes.size() ; g++)
			{
				tn = (TreeNode) sphereNodes.get(g);
				if(tn.source.isChiralCenter && !alreadyHandledCenters.contains(tn.source.atom)){
					IAtom firstInConnectedAtoms=tn.atom;
					IAtom chiralCenter = tn.source.atom;
					IAtom parent=null;
					if(tn.source.source==null)
						parent = tn.atom;
					else
						parent = tn.source.source.atom;
					
					IAtom[] sorted = null;
					List chiralNeighbours = atomContainer.getConnectedAtomsList(chiralCenter);
					if (AtomUtils.isTetrahedral(atomContainer, chiralCenter,false) > 0)
					{
						sorted = new IAtom[3];
					}
					if (AtomUtils.isTetrahedral(atomContainer, chiralCenter,false) == 1)
					{
						if (atomContainer.getBond(parent, chiralCenter).getStereo() == IBond.Stereo.DOWN || atomContainer.getBond(parent, chiralCenter).getStereo() == IBond.Stereo.DOWN_INVERTED)
						{
							boolean first=true;
							double firstangle=0;
							int firstatom=0;
							for (int k = 0; k < chiralNeighbours.size(); k++)
							{
								if (chiralNeighbours.get(k) != parent){
									if (chiralNeighbours.get(k) != parent && (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == null ||
											 atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.NONE))
									{
										double angle = BondTools.giveAngle(((IAtom) chiralNeighbours.get(k)), parent, chiralCenter);
										if(!first){
											if(angle<firstangle){
												sorted[1]=(IAtom) chiralNeighbours.get(k);
												sorted[2]=(IAtom) chiralNeighbours.get(firstatom);
											}else{
												sorted[2]=(IAtom) chiralNeighbours.get(k);
												sorted[1]=(IAtom) chiralNeighbours.get(firstatom);
											}
										}else{
											firstangle=angle;
											first=false;
											firstatom=k;
										}										
									}
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.UP)// && !isBondBroken((IAtom) chiralNeighbours.get(k), chiralCenter))
									{
										sorted[0] = (IAtom) chiralNeighbours.get(k);
									}
								}
							}
						}
						if (atomContainer.getBond(parent, chiralCenter).getStereo() == IBond.Stereo.UP || atomContainer.getBond(parent, chiralCenter).getStereo() == IBond.Stereo.UP_INVERTED)
						{
							boolean first=true;
							double firstangle=0;
							int firstatom=0;
							for (int k = 0; k < chiralNeighbours.size(); k++)
							{
								if (chiralNeighbours.get(k) != parent)
								{
									if (chiralNeighbours.get(k) != parent && (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == null ||
											 atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.NONE))
									{
										double angle = BondTools.giveAngle(((IAtom) chiralNeighbours.get(k)), parent, chiralCenter);
										if(!first){
											if(angle<firstangle){
												sorted[2]=(IAtom) chiralNeighbours.get(k);
												sorted[1]=(IAtom) chiralNeighbours.get(firstatom);
											}else{
												sorted[1]=(IAtom) chiralNeighbours.get(k);
												sorted[2]=(IAtom) chiralNeighbours.get(firstatom);
											}
										}else{
											firstangle=angle;
											first=false;
											firstatom=k;
										}										
									}
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.DOWN)// && !isBondBroken((IAtom) chiralNeighbours.get(k), chiralCenter))
									{
										sorted[0] = (IAtom) chiralNeighbours.get(k);
									}
								}
							}
						}
						if (atomContainer.getBond(parent, chiralCenter).getStereo() == CDKConstants.UNSET || atomContainer.getBond(parent, chiralCenter).getStereo() == IBond.Stereo.NONE)
						{
							boolean normalBindingIsLeft = false;
							for (int k = 0; k < chiralNeighbours.size(); k++)
							{
								if (chiralNeighbours.get(k) != parent)
								{
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == null ||
										atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.NONE)
									{
										if (BondTools.isLeft(((IAtom) chiralNeighbours.get(k)), parent, chiralCenter))
										{
											normalBindingIsLeft = true;
											break;
										}
									}
								}
							}
							for (int k = 0; k < chiralNeighbours.size(); k++)
							{
								if (chiralNeighbours.get(k) != parent)
								{
									if (normalBindingIsLeft)
									{
										if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == null ||
											atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.NONE)
										{
											sorted[0] = (IAtom) chiralNeighbours.get(k);
										}
										if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.UP || atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.UP_INVERTED)
										{
											sorted[2] = (IAtom) chiralNeighbours.get(k);
										}
										if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.DOWN || atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.DOWN_INVERTED)
										{
											sorted[1] = (IAtom) chiralNeighbours.get(k);
										}
									} else
									{
										if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.UP || atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.UP_INVERTED)
										{
											sorted[1] = (IAtom) chiralNeighbours.get(k);
										}
										if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == null ||
											atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.NONE)
										{
											sorted[0] = (IAtom) chiralNeighbours.get(k);
										}
										if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.DOWN || atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.DOWN_INVERTED)
										{
											sorted[2] = (IAtom) chiralNeighbours.get(k);
										}
									}
								}
							}
						}
					}
					if (AtomUtils.isTetrahedral(atomContainer, chiralCenter,false) == 2)
					{
						if (atomContainer.getBond(parent, chiralCenter).getStereo() == IBond.Stereo.UP)
						{
							for (int k = 0; k < chiralNeighbours.size(); k++)
							{
								if (chiralNeighbours.get(k) != parent)
								{
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.DOWN && BondTools.isLeft(((IAtom) chiralNeighbours.get(k)), parent, chiralCenter))// && !isBondBroken((IAtom) chiralNeighbours.get(k), chiralCenter))
									{
										sorted[1] = (IAtom) chiralNeighbours.get(k);
									}
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.DOWN && !BondTools.isLeft(((IAtom) chiralNeighbours.get(k)), parent, chiralCenter))// && !isBondBroken((IAtom) chiralNeighbours.get(k), chiralCenter))
									{
										sorted[2] = (IAtom) chiralNeighbours.get(k);
									}
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.UP)// && !isBondBroken((IAtom) chiralNeighbours.get(k), chiralCenter))
									{
										sorted[0] = (IAtom) chiralNeighbours.get(k);
									}
								}
							}
						}
						if (atomContainer.getBond(parent, chiralCenter).getStereo() == IBond.Stereo.DOWN)
						{
							double angle1 = 0;
							double angle2 = 0;
							IAtom atom1 = null;
							IAtom atom2 = null;
							for (int k = 0; k < chiralNeighbours.size(); k++)
							{
								if (chiralNeighbours.get(k) != parent)
								{
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.UP)// && !isBondBroken((IAtom) chiralNeighbours.get(k), chiralCenter))
									{
										if (angle1 == 0)
										{
											angle1 = BondTools.giveAngle(chiralCenter, parent, (IAtom) chiralNeighbours.get(k));
											atom1 = (IAtom) chiralNeighbours.get(k);
										} else
										{
											angle2 = BondTools.giveAngle(chiralCenter, parent, (IAtom) chiralNeighbours.get(k));
											atom2 = (IAtom) chiralNeighbours.get(k);
										}
									}
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.DOWN) // && !isBondBroken((IAtom) chiralNeighbours.get(k), chiralCenter))
									{
										sorted[1] = (IAtom) chiralNeighbours.get(k);
									}
								}
							}
							if (angle1 < angle2)
							{
								sorted[0] = atom2;
								sorted[2] = atom1;
							} else
							{
								sorted[0] = atom1;
								sorted[2] = atom2;
							}
						}
					}
					if (AtomUtils.isTetrahedral(atomContainer, chiralCenter,false) == 3)
					{
						if (atomContainer.getBond(parent, chiralCenter).getStereo() == IBond.Stereo.UP)
						{
							TreeMap hm = new TreeMap();
							for (int k = 0; k < chiralNeighbours.size(); k++)
							{
								if (chiralNeighbours.get(k) != parent)// && !isBondBroken((IAtom) chiralNeighbours.get(k), chiralCenter))
								{
									hm.put(new Double(BondTools.giveAngle(chiralCenter, parent, ((IAtom) chiralNeighbours.get(k)))), Integer.valueOf(k));
								}
							}
							Object[] ohere = hm.values().toArray();
							for (int k = ohere.length - 1; k > -1;k--)
							{
								sorted[k] = ((IAtom) chiralNeighbours.get(((Integer) ohere[k]).intValue()));
							}
						}
						if (atomContainer.getBond(parent, chiralCenter).getStereo() == null ||
							atomContainer.getBond(parent, chiralCenter).getStereo() == IBond.Stereo.NONE)
						{
							double angle1 = 0;
							double angle2 = 0;
							IAtom atom1 = null;
							IAtom atom2 = null;
							for (int k = 0; k < chiralNeighbours.size(); k++)
							{
								if (chiralNeighbours.get(k) != parent)
								{
									if ((atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == null ||
										 atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.NONE))// &&	!isBondBroken((IAtom) chiralNeighbours.get(k), chiralCenter))
									{
										if (angle1 == 0)
										{
											angle1 = BondTools.giveAngle(chiralCenter, parent, (IAtom) chiralNeighbours.get(k));
											atom1 = (IAtom) chiralNeighbours.get(k);
										} else
										{
											angle2 = BondTools.giveAngle(chiralCenter, parent, (IAtom) chiralNeighbours.get(k));
											atom2 = (IAtom) chiralNeighbours.get(k);
										}
									}
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.UP)// && !isBondBroken((IAtom) chiralNeighbours.get(k), chiralCenter))
									{
										sorted[0] = (IAtom) chiralNeighbours.get(k);
									}
								}
							}
							if (angle1 < angle2)
							{
								sorted[1] = atom2;
								sorted[2] = atom1;
							} else
							{
								sorted[1] = atom1;
								sorted[2] = atom2;
							}
						}
					}
					if (AtomUtils.isTetrahedral(atomContainer, chiralCenter,false) == 4)
					{
						if (atomContainer.getBond(parent, chiralCenter).getStereo() == IBond.Stereo.DOWN)
						{
							TreeMap hm = new TreeMap();
							for (int k = 0; k < chiralNeighbours.size(); k++)
							{
								if (chiralNeighbours.get(k) != parent) // && !isBondBroken((IAtom) chiralNeighbours.get(k), chiralCenter))
								{
									hm.put(new Double(BondTools.giveAngle(chiralCenter, parent, ((IAtom) chiralNeighbours.get(k)))), Integer.valueOf(k));
								}
							}
							Object[] ohere = hm.values().toArray();
							for (int k = ohere.length - 1; k > -1; k--)
							{
								sorted[k] = ((IAtom) chiralNeighbours.get(((Integer) ohere[k]).intValue()));
							}
						}
						if (atomContainer.getBond(parent, chiralCenter).getStereo() == null ||
							atomContainer.getBond(parent, chiralCenter).getStereo() == IBond.Stereo.NONE)
						{
							double angle1 = 0;
							double angle2 = 0;
							IAtom atom1 = null;
							IAtom atom2 = null;
							for (int k = 0; k < chiralNeighbours.size(); k++)
							{
								if (chiralNeighbours.get(k) != parent)
								{
									if ((atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == null ||
										 atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.NONE)) // &&!isBondBroken((IAtom) chiralNeighbours.get(k), chiralCenter))
									{
										if (angle1 == 0)
										{
											angle1 = BondTools.giveAngle(chiralCenter, parent, (IAtom) chiralNeighbours.get(k));
											atom1 = (IAtom) chiralNeighbours.get(k);
										} else
										{
											angle2 = BondTools.giveAngle(chiralCenter, parent, (IAtom) chiralNeighbours.get(k));
											atom2 = (IAtom) chiralNeighbours.get(k);
										}
									}
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.DOWN)// && !isBondBroken((IAtom) chiralNeighbours.get(k), chiralCenter))
									{
										sorted[2] = (IAtom) chiralNeighbours.get(k);
									}
								}
							}
							if (angle1 < angle2)
							{
								sorted[1] = atom2;
								sorted[0] = atom1;
							} else
							{
								sorted[1] = atom1;
								sorted[0] = atom2;
							}
						}
					}
					if (AtomUtils.isTetrahedral(atomContainer, chiralCenter,false) == 5)
					{
						if (atomContainer.getBond(parent, chiralCenter).getStereo() == IBond.Stereo.DOWN)
						{
							for (int k = 0; k < chiralNeighbours.size(); k++)
							{
								if (chiralNeighbours.get(k) != parent)
								{
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.UP)
									{
										sorted[0] = (IAtom) chiralNeighbours.get(k);
									}
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == null ||
										atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.NONE)
									{
										sorted[2] = (IAtom) chiralNeighbours.get(k);
									}
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.DOWN)
									{
										sorted[1] = (IAtom) chiralNeighbours.get(k);
									}
								}
							}
						}
						if (atomContainer.getBond(parent, chiralCenter).getStereo() == IBond.Stereo.UP)
						{
							for (int k = 0; k < chiralNeighbours.size(); k++)
							{
								if (chiralNeighbours.get(k) != parent)
								{
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.DOWN && BondTools.isLeft(((IAtom) chiralNeighbours.get(k)), parent, chiralCenter))// && !isBondBroken((IAtom) chiralNeighbours.get(k), chiralCenter))
									{
										sorted[0] = (IAtom) chiralNeighbours.get(k);
									}
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.DOWN && !BondTools.isLeft(((IAtom) chiralNeighbours.get(k)), parent, chiralCenter))// && !isBondBroken((IAtom) chiralNeighbours.get(k), chiralCenter))
									{
										sorted[2] = (IAtom) chiralNeighbours.get(k);
									}
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == null ||
										atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.NONE)
									{
										sorted[1] = (IAtom) chiralNeighbours.get(k);
									}
								}
							}
						}
						if (atomContainer.getBond(parent, chiralCenter).getStereo() == CDKConstants.UNSET || atomContainer.getBond(parent, chiralCenter).getStereo() == IBond.Stereo.NONE)
						{
							for (int k = 0; k < chiralNeighbours.size(); k++)
							{
								if (chiralNeighbours.get(k) != parent)
								{
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.DOWN && BondTools.isLeft(((IAtom) chiralNeighbours.get(k)), parent, chiralCenter))// && !isBondBroken((IAtom) chiralNeighbours.get(k), chiralCenter))
									{
										sorted[0] = (IAtom) chiralNeighbours.get(k);
									}
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.DOWN && !BondTools.isLeft(((IAtom) chiralNeighbours.get(k)), parent, chiralCenter))// && !isBondBroken((IAtom) chiralNeighbours.get(k), chiralCenter))
									{
										sorted[2] = (IAtom) chiralNeighbours.get(k);
									}
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.UP)
									{
										sorted[1] = (IAtom) chiralNeighbours.get(k);
									}
								}
							}
						}
					}
					if (AtomUtils.isTetrahedral(atomContainer, chiralCenter,false) == 6)
					{
						if (atomContainer.getBond(parent, chiralCenter).getStereo() == IBond.Stereo.UP)
						{
							for (int k = 0; k < chiralNeighbours.size(); k++)
							{
								if (chiralNeighbours.get(k) != parent)
								{
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.UP)
									{
										sorted[0] = (IAtom) chiralNeighbours.get(k);
									}
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == null ||
										atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.NONE)
									{
										sorted[2] = (IAtom) chiralNeighbours.get(k);
									}
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.DOWN)
									{
										sorted[1] = (IAtom) chiralNeighbours.get(k);
									}
								}
							}
						}
						if (atomContainer.getBond(parent, chiralCenter).getStereo() == IBond.Stereo.DOWN)
						{
							for (int k = 0; k < chiralNeighbours.size(); k++)
							{
								if (chiralNeighbours.get(k) != parent)
								{
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.UP && BondTools.isLeft(((IAtom) chiralNeighbours.get(k)), parent, chiralCenter))// && !isBondBroken((IAtom) chiralNeighbours.get(k), chiralCenter))
									{
										sorted[2] = (IAtom) chiralNeighbours.get(k);
									}
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.UP && !BondTools.isLeft(((IAtom) chiralNeighbours.get(k)), parent, chiralCenter))// && !isBondBroken((IAtom) chiralNeighbours.get(k), chiralCenter))
									{
										sorted[0] = (IAtom) chiralNeighbours.get(k);
									}
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == null ||
										atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.NONE)
									{
										sorted[1] = (IAtom) chiralNeighbours.get(k);
									}
								}
							}
						}
						if (atomContainer.getBond(parent, chiralCenter).getStereo() == CDKConstants.UNSET || atomContainer.getBond(parent, chiralCenter).getStereo() == IBond.Stereo.NONE)
						{
							for (int k = 0; k < chiralNeighbours.size(); k++)
							{
								if (chiralNeighbours.get(k) != parent)
								{
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.UP && BondTools.isLeft(((IAtom) chiralNeighbours.get(k)), parent, chiralCenter))// && !isBondBroken((IAtom) chiralNeighbours.get(k), chiralCenter))
									{
										sorted[2] = (IAtom) chiralNeighbours.get(k);
									}
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.UP && !BondTools.isLeft(((IAtom) chiralNeighbours.get(k)), parent, chiralCenter))// && !isBondBroken((IAtom) chiralNeighbours.get(k), chiralCenter))
									{
										sorted[0] = (IAtom) chiralNeighbours.get(k);
									}
									if (atomContainer.getBond((IAtom) chiralNeighbours.get(k), chiralCenter).getStereo() == IBond.Stereo.DOWN)
									{
										sorted[1] = (IAtom) chiralNeighbours.get(k);
									}
								}
							}
						}
					}
					if (BondTools.isSquarePlanar(atomContainer, chiralCenter))
					{
						sorted = new IAtom[3];
						//This produces a U=SP1 order in every case
						TreeMap hm = new TreeMap();
						for (int k = 0; k < chiralNeighbours.size(); k++)
						{
							if (chiralNeighbours.get(k) != parent)// && !isBondBroken((IAtom) chiralNeighbours.get(k), chiralCenter))
							{
								hm.put(new Double(BondTools.giveAngle(chiralCenter, parent, ((IAtom) chiralNeighbours.get(k)))), Integer.valueOf(k));
							}
						}
						Object[] ohere = hm.values().toArray();
						for (int k = 0; k < ohere.length; k++)
						{
							sorted[k] = ((IAtom) chiralNeighbours.get(((Integer) ohere[k]).intValue()));
						}
					}
					if (BondTools.isTrigonalBipyramidalOrOctahedral(atomContainer, chiralCenter)!=0)
					{
						sorted = new IAtom[atomContainer.getConnectedAtomsCount(chiralCenter) - 1];
						TreeMap hm = new TreeMap();
						if (atomContainer.getBond(parent, chiralCenter).getStereo() == IBond.Stereo.UP)
						{
							for (int k = 0; k < chiralNeighbours.size(); k++)
							{
								if (atomContainer.getBond(chiralCenter, (IAtom) chiralNeighbours.get(k)).getStereo() == null ||
									atomContainer.getBond(chiralCenter, (IAtom) chiralNeighbours.get(k)).getStereo() == IBond.Stereo.NONE)
								{
									hm.put(new Double(BondTools.giveAngle(chiralCenter, parent, ((IAtom) chiralNeighbours.get(k)))), Integer.valueOf(k));
								}
								if (atomContainer.getBond(chiralCenter, (IAtom) chiralNeighbours.get(k)).getStereo() == IBond.Stereo.DOWN)
								{
									sorted[sorted.length - 1] = (IAtom) chiralNeighbours.get(k);
								}
							}
							Object[] ohere = hm.values().toArray();
							for (int k = 0; k < ohere.length; k++)
							{
								sorted[k] = ((IAtom) chiralNeighbours.get(((Integer) ohere[k]).intValue()));
							}
						}
						if (atomContainer.getBond(parent, chiralCenter).getStereo() == IBond.Stereo.DOWN)
						{
							for (int k = 0; k < chiralNeighbours.size(); k++)
							{
								if (atomContainer.getBond(chiralCenter, (IAtom) chiralNeighbours.get(k)).getStereo() == null ||
									atomContainer.getBond(chiralCenter, (IAtom) chiralNeighbours.get(k)).getStereo() == IBond.Stereo.NONE)
								{
									hm.put(new Double(BondTools.giveAngle(chiralCenter, parent, ((IAtom) chiralNeighbours.get(k)))), Integer.valueOf(k));
								}
								if (atomContainer.getBond(chiralCenter, (IAtom) chiralNeighbours.get(k)).getStereo() == IBond.Stereo.UP)
								{
									sorted[sorted.length - 1] = (IAtom) chiralNeighbours.get(k);
								}
							}
							Object[] ohere = hm.values().toArray();
							for (int k = 0; k < ohere.length; k++)
							{
								sorted[k] = ((IAtom) chiralNeighbours.get(((Integer) ohere[k]).intValue()));
							}
						}
						if (atomContainer.getBond(parent, chiralCenter).getStereo() == null ||
							atomContainer.getBond(parent, chiralCenter).getStereo() == IBond.Stereo.NONE)
						{
							for (int k = 0; k < chiralNeighbours.size(); k++)
							{
								if (chiralNeighbours.get(k) != parent)
								{
									if (atomContainer.getBond(chiralCenter, (IAtom) chiralNeighbours.get(k)).getStereo() == null ||
										atomContainer.getBond(chiralCenter, (IAtom) chiralNeighbours.get(k)).getStereo() == IBond.Stereo.NONE)
									{
										hm.put(new Double((BondTools.giveAngleFromMiddle(chiralCenter, parent, ((IAtom) chiralNeighbours.get(k))))), Integer.valueOf(k));
									}
									if (atomContainer.getBond(chiralCenter, (IAtom) chiralNeighbours.get(k)).getStereo() == IBond.Stereo.UP)
									{
										sorted[0] = (IAtom) chiralNeighbours.get(k);
									}
									if (atomContainer.getBond(chiralCenter, (IAtom) chiralNeighbours.get(k)).getStereo() == IBond.Stereo.DOWN)
									{
										sorted[sorted.length - 2] = (IAtom) chiralNeighbours.get(k);
									}
								}
							}
							Object[] ohere = hm.values().toArray();
							sorted[sorted.length - 1] = ((IAtom) chiralNeighbours.get(((Integer) ohere[ohere.length - 1]).intValue()));
							if (ohere.length == 2)
							{
								sorted[sorted.length - 3] = ((IAtom) chiralNeighbours.get(((Integer) ohere[0]).intValue()));
								if (BondTools.giveAngleFromMiddle(chiralCenter, parent, ((IAtom) chiralNeighbours.get(((Integer) ohere[1]).intValue()))) < 0)
								{
									IAtom dummy = sorted[sorted.length - 2];
									sorted[sorted.length - 2] = sorted[0];
									sorted[0] = dummy;
								}
							}
							if (ohere.length == 3)
							{
								sorted[sorted.length - 3] = sorted[sorted.length - 2];
								sorted[sorted.length - 2] = ((IAtom) chiralNeighbours.get(((Integer) ohere[ohere.length - 2]).intValue()));
								sorted[sorted.length - 4] = ((IAtom) chiralNeighbours.get(((Integer) ohere[ohere.length - 3]).intValue()));
							}
						}
					}
					//This builds an onew[] containing the objects after the center of the chirality in the order given by sorted[]
					if (sorted != null)
					{
						int numberOfAtoms = 3;
						if (BondTools.isTrigonalBipyramidalOrOctahedral(atomContainer, chiralCenter)!=0)
						{
							numberOfAtoms = atomContainer.getConnectedAtomsCount(chiralCenter) - 1;
						}
						/*Object[] omy = new Object[numberOfAtoms];
						Object[] onew = new Object[numberOfAtoms];
						for (int k = getRingOpenings(atom, null).size(); k < numberOfAtoms; k++)
						{
							if (positionInVector + 1 + k - getRingOpenings(atom, null).size() < v.size())
							{
								omy[k] = v.get(positionInVector + 1 + k - getRingOpenings(atom, null).size());
							}
						}
						for (int k = 0; k < sorted.length; k++)
						{
							if (sorted[k] != null)
							{
								for (int m = 0; m < omy.length; m++)
								{
									if (omy[m] instanceof IAtom)
									{
										if (omy[m] == sorted[k])
										{
											onew[k] = omy[m];
										}
									} else
									{
										if (omy[m] == null)
										{
											onew[k] = null;
										} else
										{
											if (((List) omy[m]).get(0) == sorted[k])
											{
												onew[k] = omy[m];
											}
										}
									}
								}
							} else
							{
								onew[k] = null;
							}
						}*/
						//This is a workaround for 3624.MOL.2 I don't have a better solution currently
						/*boolean doubleentry = false;
						for (int m = 0; m < onew.length; m++)
						{
							for (int k = 0; k < onew.length; k++)
							{
								if (m != k && onew[k] == onew[m])
								{
									doubleentry = true;
								}
							}
						}
						if (!doubleentry)
						{*/
							//Make sure that the first atom in onew is the first one in the original smiles order. This is important to have a canonical smiles.
							/*if (positionInVector + 1 < v.size())
							{
								Object atomAfterCenterInOriginalSmiles = v.get(positionInVector + 1);*/
								int l = 0;
								while (sorted[0] != firstInConnectedAtoms)
								{
									IAtom placeholder = sorted[sorted.length - 1];
									for (int k = sorted.length - 2; k > -1; k--)
									{
										sorted[k + 1] = sorted[k];
									}
									sorted[0] = placeholder;
									l++;
									if (l > sorted.length)
									{
										break;
									}
								}
								for(int k=0;k<sorted.length;k++){
									IAtom toSort=sorted[k];
									for (int m = 0;  m< sphereNodes.size() ; m++)
									{
										TreeNode treeNode = (TreeNode) sphereNodes.get(m);
										if(treeNode.atom==toSort){
											treeNode.score=k*100000;
										}
									}
								}
							//}
							//This cares about ring openings. Here the ring closure (represendted by a figure) must be the first atom. In onew the closure is null.
							/*if (getRingOpenings(atom, null).size() > 0)
							{
								int l = 0;
								while (onew[0] != null)
								{
									Object placeholder = onew[0];
									for (int k = 1; k < onew.length; k++)
									{
										onew[k - 1] = onew[k];
									}
									onew[onew.length - 1] = placeholder;
									l++;
									if (l > onew.length)
									{
										break;
									}
								}
							}
							//The last in onew is a vector: This means we need to exchange the rest of the original smiles with the rest of this vector.
							if (onew[numberOfAtoms - 1] instanceof List)
							{
								for (int i = 0; i < numberOfAtoms; i++)
								{
									if (onew[i] instanceof IAtom)
									{
										List vtemp = new Vector();
										vtemp.add(onew[i]);
										for (int k = positionInVector + 1 + numberOfAtoms; k < v.size(); k++)
										{
											vtemp.add(v.get(k));
										}
										onew[i] = vtemp;
										for (int k = v.size() - 1; k > positionInVector + 1 + numberOfAtoms - 1; k--)
										{
											v.remove(k);
										}
										for (int k = 1; k < ((List) onew[numberOfAtoms - 1]).size(); k++)
										{
											v.add(((List) onew[numberOfAtoms - 1]).get(k));
										}
										onew[numberOfAtoms - 1] = ((List) onew[numberOfAtoms - 1]).get(0);
										break;
									}
								}
							}
							//Put the onew objects in the original Vector
							int k = 0;
							for (int m = 0; m < onew.length; m++)
							{
								if (onew[m] != null)
								{
									v.set(positionInVector + 1 + k, onew[m]);
									k++;
								}
							}
						}*/
					}					
				}
			}
			sortNodesByScore(sphereNodes);
		}
		for (int f = 0; f < maxSphere; f++)
		{
			sphereNodes = spheres[f];
			for (int g = 0; g < sphereNodes.size() ; g++)
			{
				tn = (TreeNode) sphereNodes.get(g);
        String localscore=tn.score+"";
        while(localscore.length()<6){
          localscore="0"+localscore;
        }
        tn.stringscore=tn.source.stringscore+""+localscore;
			}
			sortNodesByScore(sphereNodes);
		}
		HOSECode.append(centerCode);
		for (int f = 0; f < maxSphere; f++)
		{
			sphere = f + 1;
			sphereNodes = spheres[f];
			String s = getSphereCode(sphereNodes);
			HOSECode.append(s);
		}
	}

	/**
	 *  Generates the string code for a given sphere.
	 *
	 *@param  sphereNodes A vector of TreeNodes for which a string code is to be generated
	 *@return The SphereCode value
	 *@exception  org.openscience.cdk.exception.CDKException  Thrown if something goes wrong
	 */
	private String getSphereCode(List<TreeNode> sphereNodes) throws CDKException
	{
		if (sphereNodes == null || sphereNodes.size() < 1)
		{
			return sphereDelimiters[sphere - 1];
		}
		TreeNode treeNode = null;
		StringBuffer code = new StringBuffer();
		/*
		 *  append the tree node code to the HOSECode in
		 *  their now determined order, using commas to
		 *  separate nodes from different branches
		 */
		IAtom branch = sphereNodes.get(0).source.atom;
		StringBuffer tempCode = null;
		for (int i = 0; i < sphereNodes.size(); i++)
		{
			treeNode = sphereNodes.get(i);
			tempCode = new StringBuffer();
			boolean isChiral=false;
			if(i==0 && i<sphereNodes.size()-1 && sphereNodes.get(0).source.isChiralCenter){
				code.append("@");
				isChiral=true;
			}
			if (!treeNode.source.stopper && treeNode.source.atom != branch)
			{
				branch = treeNode.source.atom;
				code.append(",");
				if(i<sphereNodes.size()-1 && sphereNodes.get(i).source.isChiralCenter){
					code.append("@");
					isChiral=true;
				}
			}
			
			if (!treeNode.source.stopper && treeNode.source.atom == branch)
			{
				if (treeNode.bondType <= 4)
				{
					tempCode.append(bondSymbols[(int) treeNode.bondType]);
				} else
				{
					throw new CDKException("Unknown bond type");
				}
				if(treeNode.atom!=null && treeNode.source!=null && treeNode.source.source!=null 
						&& treeNode.source.source.source!=null
						&& (treeNode.source.atom.getSymbol().equals("C") || treeNode.source.atom.getSymbol().equals("N"))
						&& treeNode.source.bondType==2
						&& (treeNode.source.source.atom.getSymbol().equals("C") || treeNode.source.source.atom.getSymbol().equals("N"))){
					double angle1=Math.abs(BondTools.giveAngleBothMethods(treeNode.source.atom, treeNode.atom, treeNode.source.source.atom, false));
					if(angle1>Math.PI)
						angle1=2*Math.PI-angle1;
					double angle2=Math.abs(BondTools.giveAngleBothMethods(treeNode.source.source.atom, treeNode.source.atom, treeNode.source.source.source.atom,false));
					if(angle2>Math.PI)
						angle2=2*Math.PI-angle2;
					if(angle1<Math.PI*.9
							&& angle2<Math.PI*.9){
						if(isCisTrans(treeNode.atom, treeNode.source.atom,treeNode.source.source.atom, treeNode.source.source.source.atom, atomContainer)){
							tempCode.append("|");
						}else{
							tempCode.append("\\");
						}
					}
				}
				if (treeNode.atom != null && !treeNode.atom.getFlag(CDKConstants.VISITED))
				{
					tempCode.append(getElementSymbol(treeNode.symbol));
				}
				else if (treeNode.atom != null && treeNode.atom.getFlag(CDKConstants.VISITED))
				{
					tempCode.append("&");
					treeNode.stopper = true;
				}
        code.append(tempCode+createChargeCode(treeNode.atom));
				treeNode.hSymbol = tempCode.toString();
			}
			if (treeNode.atom != null) treeNode.atom.setFlag(CDKConstants.VISITED, true);
			if (treeNode.source.stopper) treeNode.stopper = true;
		}
		code.append(sphereDelimiters[sphere - 1]);
		return code.toString();
	}


	/**
	 *  Gets the element rank for a given element symbol as given in Bremser's
	 *  publication.
	 *
	 *@param  symbol  The element symbol for which the rank is to be determined
	 *@return         The element rank
	 */
	private double getElementRank(String symbol)
	{
		for (int f = 0; f < rankedSymbols.length; f++)
		{
			if (rankedSymbols[f].equals(symbol))
			{
				return symbolRankings[f];
			}
		}
		IIsotope isotope = isotopeFac.getMajorIsotope(symbol);
		if (isotope.getMassNumber() != null) {
			return ((double) 800000 - isotope.getMassNumber());
		}
        return 800000;
	}

	/**
	 *  Returns the Bremser-compatible symbols for a given element. Silicon, for
	 *  example, is actually "Q". :-)
	 *
	 *@param  sym  The element symbol to be converted
	 *@return      The converted symbol
	 */
	private String getElementSymbol(String sym)
	{
		if (sym.equals("Si"))
		{
			return "Q";
		}
		if (sym.equals("Cl"))
		{
			return "X";
		}
		if (sym.equals("Br"))
		{
			return "Y";
		}
		if (sym.equals(","))
		{
			return "";
		}
		return sym;
	}


	/**
	 *  Determines the ranking score for each node, allowing for a sorting of nodes
	 *  within one sphere.
	 *
	 *@param  sphereNodes The nodes for which the score is to be calculated.
	 *@exception  org.openscience.cdk.exception.CDKException  Thrown if something goes wrong.
	 */
	private void calculateNodeScores(List<TreeNode> sphereNodes) throws CDKException
	{
		TreeNode treeNode = null;
		for (int i = 0; i < sphereNodes.size(); i++)
		{
			treeNode = (TreeNode) sphereNodes.get(i);
			treeNode.score += getElementRank(treeNode.symbol);
			if (treeNode.bondType <= 4)
			{
				treeNode.score += bondRankings[(int) treeNode.bondType];
			} else
			{
				throw new CDKException("Unknown bond type encountered in HOSECodeGenerator");
			}
		}
	}

	/**
	 *  Sorts the nodes (atoms) in the sphereNode vector according to their score.
	 *  This is used for the essential ranking of nodes in HOSE code sphere.
	 *
	 *@param  sphereNodes  A vector with sphere nodes to be sorted.
	 */
	private void sortNodesByScore(List<TreeNode> sphereNodes)
	{
		TreeNode obj;
		boolean changed;
		if (sphereNodes.size() == 0) return;
		/*
		 *  Now we sort by score
		 */
		do
		{
			changed = false;
			for (int i = 0; i < sphereNodes.size() - 1; i++)
			{
				if (((TreeNode) sphereNodes.get(i + 1)).stringscore.compareTo(((TreeNode) sphereNodes.get(i)).stringscore)>0)
				{
					obj = sphereNodes.get(i + 1);
					sphereNodes.remove(i + 1);
					sphereNodes.add(i, obj);
					changed = true;
				}
			}
		} while (changed);
		/* Having sorted a sphere, we label the nodes with their sort order */
		TreeNode temp = null;
		for (int i = 0; i < sphereNodes.size(); i++)
		{
			temp = ((TreeNode) sphereNodes.get(i));
			temp.sortOrder = sphereNodes.size() - i;
		}
	}

	/**
	 *  If we use less than four sphere, this fills up the code with the missing
	 *  delimiters such that we are compatible with Bremser's HOSE code table.
	 */
	private void fillUpSphereDelimiters()
	{
		logger.debug("Sphere: " + sphere);
		for (int f = sphere; f < 4; f++)
		{
			HOSECode.append(sphereDelimiters[f]);
		}
	}

	
	class TreeNodeComparator implements Comparator<TreeNode> {
    /**
     *The compare method, compares by canonical label of atoms
     *
     * @param  obj1  The first TreeNode
     * @param  obj2  The second TreeNode
     * @return       -1,0,1
     */
    public int compare(TreeNode obj1, TreeNode obj2) {
    	if(obj1==null || obj2==null || ((TreeNode) obj1).getAtom()==null || ((TreeNode) obj2).getAtom()==null)
    		return 0;
    	Long label1 = (Long)((TreeNode) obj1).getAtom().getProperty(InvPair.CANONICAL_LABEL);
    	Long label2 = (Long)((TreeNode) obj2).getAtom().getProperty(InvPair.CANONICAL_LABEL);
    	if(label1==null || label2==null)
    		return 0;
    	if (label1.intValue() < label2.intValue()) {
    		return (-1);
    	}
    	if (label1.intValue() > label2.intValue()) {
    		return (1);
    	}
    	return (0);
    }
  }
    
  /**
	 *  Helper class for storing the properties of a node in our breadth first
	 *  search.
	 *
	 * @author     steinbeck
	 * @cdk.created    2002-11-16
	 */
	class TreeNode
	{
		String symbol;
		TreeNode source;
		IAtom atom;
		double bondType;
		int degree;
		long score;
		int ranking;
		int sortOrder = 1;
		List<TreeNode> childs = null;
		String hSymbol = null;
		boolean stopper = false;
    String stringscore="";
    IBond.Stereo stereo;
    boolean isChiralCenter=false;

		/**
		 *  Constructor for the TreeNode object.
		 *
		 *@param  symbol    The Element symbol of the node
		 *@param  source    The preceding node for this node
		 *@param  atom      The IAtom object belonging to this node
		 *@param  bondType  The bond type by which this node was connect to its
		 *      predecessor
		 *@param  score     The score used to rank this node within its sphere.
		 *@param  degree    Description of the Parameter
		 */
		TreeNode(String symbol, TreeNode source, IAtom atom, double bondType, int degree, long score, IBond.Stereo stereo, boolean isChiralCenter)
		{
			this.symbol = symbol;
			this.source = source;
			this.atom = atom;
			this.degree = degree;
			this.score = score;
			this.bondType = bondType;
			this.stereo = stereo;
			this.isChiralCenter = isChiralCenter;
			ranking = 0;
			sortOrder = 1;
			childs = new ArrayList<TreeNode>();
		}
    
    public IAtom getAtom(){
      return atom;
    }


		/**
		 *  A TreeNode is equal to another TreeNode if it
		 *  stands for the same atom object.
		 *
		 *@param  o  The object that we compare this TreeNode to
		 *@return    True, if the this TreeNode's atom object equals the one of the other TreeNode
		 */
		public boolean equals(Object o)
		{
			try
			{
				if (this.atom == ((TreeNode) o).atom)
				{
					return true;
				}
			}
			catch(Exception exc)
			{
				/* we do nothing here because anything 
				that could seriously happen here is the we 
				got something which is not a TreeNode and then got 
				a class cast exception. Thus we can just wait until the
				end of the method returns a "false" */
			}
			return false;
		}
		
		public String toString()
		{
			String s = "";
			try
			{
				s += (atomContainer.getAtomNumber(atom) + 1);
				s += " " + hSymbol;
				s += "; s=" + score; 
				s += "; r=" + ranking;
				s += "; d = " + degree;
			}
			catch(Exception exc)
			{
				return exc.toString();	
			}
			return s;
		}
	}
  
  
	public List<IAtom> getNodesInSphere(int sphereNumber){
		sphereNodes = spheres[sphereNumber-1];
		List<IAtom> atoms=new ArrayList<IAtom>();
		for (int g = 0; g < sphereNodes.size() ; g++) {
			atoms.add(((TreeNode) sphereNodes.get(g)).atom);
		}
		return(atoms);
	}
	
	  /**
	   *  Says if two atoms are in cis or trans position around a double bond.
	   *  The atoms have to be given to the method like this:  firstOuterAtom - firstInnerAtom = secondInnterAtom - secondOuterAtom
	   *
	   * @param  firstOuterAtom    See above.
	   * @param  firstInnerAtom    See above.
	   * @param  secondInnerAtom   See above.
	   * @param  secondOuterAtom   See above.
	   * @param  ac                The atom container the atoms are in.
	   * @return                   true=trans, false=cis.
	   * @exception  CDKException  The atoms are not in a double bond configuration (no double bond in the middle, same atoms on one side)
	   */
	  public static boolean isCisTrans(IAtom firstOuterAtom, IAtom firstInnerAtom, IAtom secondInnerAtom, IAtom secondOuterAtom, IAtomContainer ac) throws CDKException {
	    boolean firstDirection = BondTools.isLeft(secondOuterAtom, firstInnerAtom, secondInnerAtom);
	    boolean secondDirection = BondTools.isLeft(firstOuterAtom, secondInnerAtom, firstInnerAtom);
	      return firstDirection == secondDirection;
	  }

}


