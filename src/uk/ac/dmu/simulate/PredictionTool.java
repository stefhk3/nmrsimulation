package uk.ac.dmu.simulate;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.TreeMap;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.HOSECodeGenerator;

public class PredictionTool
{
  public Map<String, Map<String, PredictionTool.ValueBean>> mapsmap = new TreeMap<String, Map<String, PredictionTool.ValueBean>>();
  
  public Map<String, float[]> usedHoseCodes = new TreeMap<String, float[]>();
  




  public PredictionTool()
    throws IOException
  {
    String filename = "nmrshiftdbc.csv";
    InputStream ins = getClass().getClassLoader().getResourceAsStream(filename);
    BufferedReader reader = new BufferedReader(new InputStreamReader(ins));
    String input;
    while ((input = reader.readLine()) != null) { 
    	if(!input.equals("///"))
      if (!input.startsWith("#")) {
        StringTokenizer st2 = new StringTokenizer(input, "_");
        st2.nextToken();
        String solvent = st2.nextToken();
        String code = st2.nextToken();
        Double min = new Double(st2.nextToken());
        Double av = new Double(st2.nextToken());
        Double max = new Double(st2.nextToken());
        Integer numValues = new Integer(st2.nextToken());
        if(!mapsmap.containsKey(solvent))
        	mapsmap.put(solvent,  new HashMap<String,PredictionTool.ValueBean>());
        mapsmap.get(solvent).put(code, new PredictionTool.ValueBean(min.floatValue(), av.floatValue(), max.floatValue(), numValues.intValue()));
      }
    }
    reader.close();
    filename = "nmrshiftdbh.csv";
    ins = getClass().getClassLoader().getResourceAsStream(filename);
    reader = new BufferedReader(new InputStreamReader(ins));
    while ((input = reader.readLine()) != null) { 
    	if(!input.equals("///"))
      if (!input.startsWith("#")) {
        StringTokenizer st2 = new StringTokenizer(input, "_");
        st2.nextToken();
        String solvent = st2.nextToken();
        String code = st2.nextToken();
        Double min = new Double(st2.nextToken());
        Double av = new Double(st2.nextToken());
        Double max = new Double(st2.nextToken());
        Integer numValues = new Integer(st2.nextToken());
        if(!mapsmap.containsKey(solvent))
        	mapsmap.put(solvent,  new HashMap<String,PredictionTool.ValueBean>());
        mapsmap.get(solvent).put(code, new PredictionTool.ValueBean(min.floatValue(), av.floatValue(), max.floatValue(), numValues.intValue()));
      }
    }
    reader.close();
  }
  


















































  public float[] generalPredict(IAtomContainer mol, IAtom a, boolean calculated, boolean measured, boolean hose3d, int ignoreSpectrum, int ignoreSpectrumEnd, StringBuffer comment, boolean commentWithMinMax, boolean withRange, Map<Double, double[]> predictionValuesForApplet, int maxSpheresToUse, boolean cache, StringBuffer hoseCodeOut, int spheresMax, boolean fromDB, boolean trueonly, String solvent)
    throws CDKException, Exception
  {
    HOSECodeGenerator hcg = new HOSECodeGenerator();
    if ((maxSpheresToUse == -1) || (maxSpheresToUse > spheresMax)) {
      maxSpheresToUse = spheresMax;
    }
    float[] returnValues = new float[9];
    
    String hoseCodeOrigExtended = null;
    String hoseCode = null;
    synchronized (mol) {
      hoseCode = hcg.getHOSECode(mol, a, maxSpheresToUse, false);
      if (hose3d) {
        ExtendedHOSECodeGenerator ehcg = new ExtendedHOSECodeGenerator();
        hoseCodeOrigExtended = ehcg.getHOSECode(mol, a, maxSpheresToUse, true);
      } else {
        hoseCodeOrigExtended = hoseCode;
      }
    }
    if ((hose3d) && (usedHoseCodes.containsKey(hoseCodeOrigExtended)))
      return (float[])usedHoseCodes.get(hoseCodeOrigExtended);
    if ((!hose3d) && (usedHoseCodes.containsKey(hoseCode))) {
      return (float[])usedHoseCodes.get(hoseCode);
    }
    for (int spheres = maxSpheresToUse; spheres > 0; spheres--) {
    	//if we are down to 3 spheres with a solvent, we better go for unreported.
    	if(spheres==3 && !solvent.equals("Unreported"))
    		return generalPredict(mol,a,calculated,measured,hose3d,ignoreSpectrum,ignoreSpectrumEnd,comment,commentWithMinMax,withRange,predictionValuesForApplet,maxSpheresToUse,cache,hoseCodeOut,spheresMax,fromDB,trueonly,"Unreported");
      StringBuffer hoseCodeBuffer = new StringBuffer();
      StringTokenizer st = new StringTokenizer(hoseCode, "()/");
      StringBuffer hoseCodeBufferExtended = null;
      StringTokenizer stExtended = null;
      if (hose3d) {
        hoseCodeBufferExtended = new StringBuffer();
        stExtended = new StringTokenizer(hoseCodeOrigExtended, "()/");
      }
      for (int k = 0; k < spheres; k++) {
        if (st.hasMoreTokens()) {
          String partcode = st.nextToken();
          hoseCodeBuffer.append(partcode);
          if (hose3d) {
            partcode = stExtended.nextToken();
            hoseCodeBufferExtended.append(partcode);
          }
        }
        if (k == 0) {
          hoseCodeBuffer.append("(");
          if (hose3d)
            hoseCodeBufferExtended.append("(");
        } else if (k == 3) {
          hoseCodeBuffer.append(")");
          if (hose3d)
            hoseCodeBufferExtended.append(")");
        } else {
          hoseCodeBuffer.append("/");
          if (hose3d) {
            hoseCodeBufferExtended.append("/");
          }
        }
      }
      PredictionTool.ValueBean l = null;
      if (hose3d)
        l = mapsmap.get(solvent).get(hoseCodeBufferExtended.toString());
    	if(solvent.equals("Unreported") && l==null)
      		l = mapsmap.get("Methanol-D4 (CD3OD)").get(hoseCodeBufferExtended.toString());
      	else if(solvent.equals("Unreported") && l==null)
      		l = mapsmap.get("Chloroform-D1 (CDCl3)").get(hoseCodeBufferExtended.toString());
      if (l != null) {
        returnValues[0] = l.min;
        returnValues[1] = l.average;
        returnValues[2] = l.max;
        returnValues[4] = spheres;
        usedHoseCodes.put(hoseCodeOrigExtended, returnValues);
        return returnValues;
      }
      l = mapsmap.get(solvent).get(hoseCodeBuffer.toString());
    	if(solvent.equals("Unreported") && l==null)
      		l = mapsmap.get("Methanol-D4 (CD3OD)").get(hoseCodeBuffer.toString());
      	else if(solvent.equals("Unreported") && l==null)
      		l = mapsmap.get("Chloroform-D1 (CDCl3)").get(hoseCodeBuffer.toString());
      if (l != null)
      {
        returnValues[0] = l.min;
        returnValues[1] = l.average;
        returnValues[2] = l.max;
        returnValues[4] = spheres;
        usedHoseCodes.put(hoseCodeOrigExtended, returnValues);
        return returnValues;
      }
    }
    
    returnValues[0] = -1.0F;
    returnValues[1] = -1.0F;
    returnValues[2] = -1.0F;
    returnValues[4] = -1.0F;
    usedHoseCodes.put(hoseCodeOrigExtended, returnValues);
    
    return returnValues;
  }
  














  public float[] predict(IAtomContainer mol, IAtom atom, boolean useExtendedHOSECode, String solvent)
    throws CDKException, Exception
  {
    return generalPredict(mol, atom, true, true, useExtendedHOSECode, -1, -1, new StringBuffer(), false, true, null, 6, false, null, 6, false, true, solvent);
  }
  
  /*public Map<Integer, Map<Integer, float[]>> predict(File sdf, boolean useExtendedHOSECode, int nrOfProcessors) throws IOException, CDKException, ClassNotFoundException {
    ExecutorService eservice = Executors.newFixedThreadPool(nrOfProcessors);
    CompletionService<Map<Integer, Map<Integer, float[]>>> cservice = new ExecutorCompletionService(eservice);
    IteratingMDLReader iter = new IteratingMDLReader(new FileInputStream(sdf), DefaultChemObjectBuilder.getInstance());
    int molcount = 0;
    while (iter.hasNext()) {
      IAtomContainer mol = iter.next();
      cservice.submit(new PredictionTool.Task(molcount, mol, this, useExtendedHOSECode));
      molcount++;
    }
    Map<Integer, Map<Integer, float[]>> result = new HashMap();
    for (int index = 0; index < molcount; index++) {
      try {
        Map<Integer, Map<Integer, float[]>> taskResult = (Map)cservice.take().get();
        int actualIndex = ((Integer)taskResult.keySet().iterator().next()).intValue();
        result.put(Integer.valueOf(actualIndex), (Map)taskResult.get(Integer.valueOf(actualIndex)));
      }
      catch (InterruptedException localInterruptedException) {}catch (ExecutionException localExecutionException) {}
    }
    
    return result;
  }*/
  

public static class ValueBean {
    public float min;
    public float average;
    public float max;
    public int numValues;
    public int numValues2;

    public ValueBean(float a, float b, float c, int numValues) {
        this.min = a;
        this.average = b;
        this.max = c;
        this.numValues = numValues;
    }

    public ValueBean(float a, float b, float c, int numValues, int numValues2) {
        this(a, b, c, numValues);
        this.numValues2 = numValues2;
    }
}
}
