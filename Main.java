package com.company;

import java.awt.image.AreaAveragingScaleFilter;
import java.math.BigInteger;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.ArrayList;

import java.util.BitSet;


public class Main {
    private static int Bfactor = 15;
    private static int errorMargin = 1;

    private static boolean[] markedRows;

    public static void main(String[] args) {


        //System.out.println(myPower(new BigInteger("" + 5), new BigInteger("" + 3)).toString());


        //BigInteger[] bla = startOffsets(BigInteger.valueOf(1232222), BigInteger.valueOf(13));
        //System.out.println(bla[0].toString() + ":" + bla[1].toString());
        //System.out.println(generateFactorBase(1232222).toString());
        //BigInteger tS =  tonelliShanks(BigInteger.valueOf(1232222), BigInteger.valueOf(13));
        //System.out.println(tS.toString());

        /*for (int i = 10000000; i < 20000000; i++) {
            BigInteger[] res = primePower(BigInteger.valueOf(i));
            if (res != null) {
                System.out.println(i + " -> " + res[0] + "^" + res[1]);
            }
            else {

               ArrayList<Integer> asd = generateFactorBase(i);

                //System.out.print(asd.toString() + " // ");
                if (asd.size() < 2) {
                    //System.out.println(i + " -> Only 2");
                }
                else {
                    BigInteger[] wall = startOffsets(BigInteger.valueOf(i), BigInteger.valueOf(asd.get(1)));
                    System.out.println(i + " -> " + wall[0] + ":" + wall[1]);
                }
            }
        }
        */

        BigInteger N = new  BigInteger("" + 1232124);
        ArrayList<Integer> factorBase = generateFactorBase(N);
        ArrayList<BigInteger> smoothNumbers = generateSmoothNumbers(N, factorBase);

        //System.out.println(smoothNumbers);
    }


    private static BigInteger myPower(BigInteger a, BigInteger b) {
        BigInteger res = BigInteger.ONE;
        while (b.compareTo(BigInteger.ZERO) > 0) {
            res = res.multiply(a);
            b = b.subtract(BigInteger.ONE);
        }
        return res;
    }

    private static BigInteger[] primePower(BigInteger N) {
       BigInteger start = N.divide(BigInteger.valueOf(2));



       for (BigInteger power = BigInteger.valueOf(2); N.compareTo(myPower(BigInteger.valueOf(2),power)) > 0; power = power.add(BigInteger.ONE)) {

           BigInteger base = start;
           BigInteger interval = start;

           while (interval.compareTo(BigInteger.ZERO) > 0) {
               //System.out.println("Base: " + base.toString() + ", Power: " + power.toString() + ", Res: " + myPower(base, power));
               int comp = myPower(base, power).compareTo(N);
               if (comp == 0) {
                   // We have a prime power!
                   return new BigInteger[]{base, power};
               }
               else if (comp < 0) {
                   // To low!
                   interval = interval.divide(BigInteger.valueOf(2));
                   base = base.add(interval);

               } else {
                   // To high!
                   interval = interval.divide(BigInteger.valueOf(2));
                   base = base.subtract(interval);

               }

            }
       }

        return null;
    }

    // SNATCHED
    public static int legendre(BigInteger a, BigInteger p) {

    }

    // SNATCHED
    public static BigInteger tonelliShanks(BigInteger a, BigInteger bigP) {

    }

    public static BigInteger[] startOffsets(BigInteger a, BigInteger bigP) {
        BigInteger tS =  tonelliShanks(a, bigP);
        BigInteger squareRoot = bigIntSqRootCeil(a);

        BigInteger x = tS.subtract(squareRoot).mod(bigP);
        BigInteger y = bigP.subtract(tS).subtract(squareRoot).mod(bigP);

        //BigInteger x = tonelliShanks(a, bigP).subtract(bigIntSqRootCeil(a)).mod(bigP);
        //BigInteger y = bigP.subtract(tonelliShanks(a,bigP)).subtract(bigIntSqRootCeil(a)).mod(bigP);
        return new BigInteger[]{x, y};
    }

    public static BigInteger bigIntSqRootCeil(BigInteger x) throws IllegalArgumentException {
        if (x.compareTo(BigInteger.ZERO) < 0) {
            throw new IllegalArgumentException("Negative argument.");
        }
        // square roots of 0 and 1 are trivial and
        // y == 0 will cause a divide-by-zero exception
        if (x == BigInteger.ZERO || x == BigInteger.ONE) {
            return x;
        } // end if
        BigInteger two = BigInteger.valueOf(2L);
        BigInteger y;
        // starting with y = x / 2 avoids magnitude issues with x squared
        for (y = x.divide(two);
             y.compareTo(x.divide(y)) > 0;
             y = ((x.divide(y)).add(y)).divide(two));
        if (x.compareTo(y.multiply(y)) == 0) {
            return y;
        } else {
            return y.add(BigInteger.ONE);
        }
    } // end bigIntSqRootCeil

    public static boolean isPrime(int n) {
        if(n < 2) return false;
        if(n == 2 || n == 3) return true;
        if(n%2 == 0 || n%3 == 0) return false;
        long sqrtN = (long)Math.sqrt(n)+1;
        for(long i = 6L; i <= sqrtN; i += 6) {
            if(n%(i-1) == 0 || n%(i+1) == 0) return false;
        }
        return true;
    }

    private static ArrayList<Integer> generateFactorBase(BigInteger N) {
        int B = Bfactor * (int)  Math.ceil(Math.exp(Math.sqrt(bigLog(N)*Math.log(bigLog(N)))/2 ));

        //int B = 1 * (int) Math.ceil(Math.exp(Math.sqrt(bigLog(N)*Math.log(bigLog(N)))));
        System.out.println("N: " + N + ", B: " + B + "\n");

        ///////////////////////////////
        // GENERATION OF FACTOR BASE //
        ///////////////////////////////

        ArrayList<Boolean> primeSieve = new ArrayList<Boolean>();
        for (int i = 0; i < B; i++) {
            primeSieve.add(i, true);
        }
        primeSieve.set(0, false);
        primeSieve.set(1, false);

        for (int i = 2; i < B; i++) {
            if (primeSieve.get(i) == true) {
                for (int j = i*i; j < B; j+=i) { // Are we going one to low with respect to log?
                    //System.out.println(i + " : " + j);
                    primeSieve.set(j, false);
                }
            }
        }

        ArrayList<Integer> factorBase = new ArrayList<Integer>();
        factorBase.add(2);
        //factorBase.add(-1);
        for (int i = 3; i < B; i++) {
            if (primeSieve.get(i) == true) {
                if (legendre(N, BigInteger.valueOf(i)) == 1)    // <- FIX FOR EFFECTIVITY!
                    factorBase.add(i);

            }
        }

        System.out.println(factorBase.toString());
        return factorBase;
    }

    private static float bigLog(BigInteger N) {
        return (float) (N.bitLength() * Math.log(2));
    }


    private static ArrayList<BigInteger> generateSmoothNumbers(BigInteger N, ArrayList<Integer> factorBase) {
        ArrayList<BigInteger> smoothNumbers = new ArrayList<BigInteger>();
        int intervalSize = factorBase.get(factorBase.size()-1);

        BigInteger Nsqrt = bigIntSqRootCeil(N);
        BigInteger firstA =  bigIntSqRootCeil(N); // BigInteger.valueOf(0);

        Integer[] primePositionsA = new Integer[factorBase.size()-1];
        Integer[] primePositionsB = new Integer[factorBase.size()-1];

        for (int i = 1; i < factorBase.size(); i++) {
            BigInteger[] currentStartOffset = startOffsets(N, BigInteger.valueOf(factorBase.get(i)));
            primePositionsA[i-1] = currentStartOffset[0].intValue();
            primePositionsB[i-1] = currentStartOffset[1].intValue();
        }

        float[] logTable = new float[factorBase.size()];
        for (int i = 0; i < factorBase.size(); i++) {
            logTable[i] = (float) Math.log(factorBase.get(i));
        }

        float[] logQ = new float[factorBase.get(factorBase.size()-1)];


        //for (int a = firstA; smoothNumbers.size() < factorBase.size()+2; a+= intervalSize) {
        while(true) {

            //for (int a = firstA; a == firstA; a+= intervalSize) {
            BigInteger b = firstA.add(BigInteger.valueOf(intervalSize)); // Interval: [a,b]


            for (int i = 0; i < intervalSize; i++) {
                logQ[i] = bigLog((Nsqrt.add(firstA).multiply(Nsqrt.add(firstA)).subtract(N))); // May be problematic...

                //logQ[i] = (float) Math.log( Math.pow(Math.floor(Math.sqrt(N.intValue())) + firstA.intValue()+i, 2) - N.intValue() );
            }

            // Take care of two's seperatly...
            for (int i = 0; i < intervalSize; i++) {
                if (firstA.add(BigInteger.valueOf(i)).mod(BigInteger.valueOf(2)).intValue() == 0) {  // <- Very experimental!
                    logQ[i] -=  logTable[0];
                }
            }


            for (int i = 0; i < factorBase.size()-1; i++) { // i = id of factor



                while (true) {
                    logQ[primePositionsA[i]] -= logTable[i+1];
                    primePositionsA[i] += factorBase.get(i+1);

                    if (primePositionsA[i] >= intervalSize) {
                        primePositionsA[i] %= intervalSize;
                        break;
                    }
                }


                while (true) {
                    logQ[primePositionsB[i]] -= logTable[i+1];
                    primePositionsB[i] += factorBase.get(i+1);

                    if (primePositionsB[i] >= intervalSize) {
                        primePositionsB[i] %= intervalSize;
                        break;
                    }
                }


                int bla = 0;
                // Add the other root...


            }

            for (int i = 0; i < intervalSize; i++) {
                if (Math.abs(logQ[i]) < errorMargin ) {

                    //System.out.println("HAT!");

                    if (bigIsSmooth(factorBase, firstA.add(BigInteger.valueOf(i)))) {
                        System.out.println(smoothNumbers.size() + " of " + (factorBase.size() +2));
                        smoothNumbers.add(firstA.add(BigInteger.valueOf(i)));

                    }

                    /*
                    if (isSmooth(factorBase, firstA.intValue() + i))
                    {
                            //System.out.print((firstA.intValue()+i) + " -> ");
                            //System.out.println("true...");
                            smoothNumbers.add(firstA.add(BigInteger.valueOf(i)));

                        System.out.println("HIT! " + smoothNumbers.size());
                        if (smoothNumbers.size() > factorBase.size() +2) {
                            System.out.println("\nComplete!");
                            return null;
                        }


                        if (smoothNumbers.size() > factorBase.size()+1 ) {
                            ArrayList<Integer> sN = new ArrayList<Integer>();
                            ArrayList<Integer> fB = new ArrayList<Integer>();
                            for (int ii = 0; ii < smoothNumbers.size(); ii++)
                                sN.add(smoothNumbers.get(ii).intValue());
                            for (int ii = 0; ii < fB.size(); ii++)
                                fB.add(smoothNumbers.get(ii).intValue());

                            int res = generateFactor(sN, fB, N.intValue());
                            if (res != 1 && res != N.intValue()){
                                System.out.println("------------------");
                                System.out.println("OLD: " + res);
                                int x = generateFactor2(smoothNumbers, factorBase, N);
                                System.out.println("NEW: " + x + "\n\n\n");
                                return null;
                            }
                        }

                    }
                    */

                    //else
                    //    System.out.println("false...");

                }
            }
            firstA = firstA.add(BigInteger.valueOf(intervalSize));

        }

        //return null;

    }

    public static int generateFactor2(ArrayList<BigInteger> smoothNumbers, ArrayList<Integer> factorBase, BigInteger N) {
        boolean[][] bitArray = new boolean[smoothNumbers.size()][factorBase.size()];

        //System.out.println("SMOOTH: " + smoothNumbers.size() + ", FACTORBASE: " + factorBase.size());

        // Set up matrix
        for (int i = 0; i < smoothNumbers.size(); i++) {
            for (int j = 0 ; j < factorBase.size(); j++) {
                bitArray[i][j] = bigHasEvenNumberOfPrime(smoothNumbers.get(i), factorBase.get(j));
                if (bitArray[i][j])
                    System.out.print("1 ");
                else
                    System.out.print("0 ");
            }
            System.out.println();
        }

        for (int i = 0;;i++) {
            boolean[] dep = gaussWrapper(bitArray, i);
            boolean isOk = false;
            System.out.println("DEP LENGTH: " + dep.length);
            for (int j = 0; j < dep.length; j++) {
                System.out.print(dep[j] + " ");
                // Check if zero...
                if (dep[j]) {
                    isOk = true;
                    break;
                }
            }
            System.out.println("\\");
            if (!isOk) {
                return 1;
            } else {
                BigInteger res = bigCalcNum(dep, smoothNumbers, N);
            }
        }
    }

    public static int generateFactor(ArrayList<Integer> smoothNumbers, ArrayList<Integer> factorBase, int N) {
        //System.out.print("Start... ");
        ArrayList<Integer> exponentVecors = new ArrayList<Integer>(smoothNumbers.size());

        for (int i = 0; i < smoothNumbers.size(); i++) {
            int currExpV = 0;
            for (int j = 0; j < factorBase.size(); j++) {
                int bla = (int)Math.pow(smoothNumbers.get(i),2) % N;

                //System.out.println("B: " + bla +", S: " + smoothNumbers.get(i) + ", F: " + factorBased.get(j) + ", R: " + hasEvenNumberOfPrime(bla, factorBased.get(j)));
                currExpV |= (int) Math.pow(2,j) * hasEvenNumberOfPrime(bla, factorBase.get(j));
            }
            exponentVecors.add(currExpV);
        }


        // Brute force all possible linear combinations...
        int validLinComb = 0;
        for (int i = 1; i < (int)Math.pow(2, smoothNumbers.size()); i++) { // Don't run over 0... will or else we will obtain a trivial solution...
            int currSum = 0;

            for (int j = 0; j < smoothNumbers.size(); j++) {
                if ((i & (int) Math.pow(2,j)) != 0){
                    for (int k = 0; k < factorBase.size(); k++) {
                        int bla = (int)Math.pow(smoothNumbers.get(j),2) % N;
                        //System.out.print(hasEvenNumberOfPrime(bla, factorBased.get(k)) + " ");
                        currSum ^= hasEvenNumberOfPrime(bla, factorBase.get(k))*(int)Math.pow(2,k);
                    }
                }
            }
            if (currSum == 0) {
                //validLinComb = i;
                //break;

                int x = calcNum(i, smoothNumbers, N);
                if (x != 1 && x != N) {
                    //System.out.println("End!");
                    return x;
                }
            }
        }
        //System.out.println("End!");
        return 1;
    }

    public static int calcNum(int validLinComb, ArrayList<Integer> smoothNumbers, int N) {
        ArrayList<Integer> chosenSmoothNumbers = new ArrayList<Integer>();
        for (int i = 0; i < smoothNumbers.size(); i++){
            if (((int)Math.pow(2,i) & validLinComb) != 0){
                chosenSmoothNumbers.add(smoothNumbers.get(i));
            }
        }

        //System.out.println("Valid smooth numbers: " + chosenSmoothNumbers.toString());

        // Generating number from chosen smooth number
        int a = 1;
        int b = 1;
        for (int i = 0; i < chosenSmoothNumbers.size(); i++) {
            a *= chosenSmoothNumbers.get(i);
            b *= (chosenSmoothNumbers.get(i)*chosenSmoothNumbers.get(i))%N;
        }
        a %= N;
        b = (int)Math.sqrt(b);

        int c = gcdThing(a-b,N);
        //System.out.println("Answer: " + c);

        return c;
    }
    public static BigInteger bigCalcNum(boolean[] linearCombination, ArrayList<BigInteger> smoothNumbers, BigInteger N) {
        BigInteger a = BigInteger.ONE;
        BigInteger b = BigInteger.ONE;

        for (int i = 0; i < linearCombination.length; i++) {
            if (linearCombination[i]) {
                a = a.multiply(smoothNumbers.get(i));
                b = b.multiply(smoothNumbers.get(i)).multiply(smoothNumbers.get(i)).mod(N);
            }
        }
        a = a.mod(N);
        b = bigIntSqRootCeil(b);   // May need to be floor....

        return (a.subtract(b)).gcd(N);
    }

    public static Boolean bigIsSmooth(ArrayList<Integer> factorBase, BigInteger N) {
        for (int i = 0; i < factorBase.size() && N.compareTo(BigInteger.ONE) > 0; i++) {
            //int currentPrime = pB.get(i);
            //if (N % currentPrime == 0) {
            //    N = N/currentPrime; //N = N % pB.get(i);
            //}

            N = bigFactorOut(N, factorBase.get(i));
        }
        if (N.equals(BigInteger.ONE)){
            return true;
        } else {
            return false;
        }
    }

    public static Boolean isSmooth(ArrayList<Integer> factorBase, int N) {
        for (int i = 0; i < factorBase.size() && N != 1; i++) {
            //int currentPrime = pB.get(i);
            //if (N % currentPrime == 0) {
            //    N = N/currentPrime; //N = N % pB.get(i);
            //}

            N = FactorOut(N, factorBase.get(i));
        }
        if (N == 1) {
            return true;
        }
        else {
            return false;
        }
    }

    public static int FactorOut(int N, int p) {
        while (N % p == 0) {
            N = N/p; //N = N % p;
        }
        return N;
    }

    public static BigInteger bigFactorOut(BigInteger N, int p) {
        while (N.mod(BigInteger.valueOf(p)).compareTo(BigInteger.ZERO) == 0) {
            N = N.divide(BigInteger.valueOf(p));
        }
        return N;
    }

    public static boolean bigHasEvenNumberOfPrime(BigInteger N, Integer p) {
        int i = 0;
        BigInteger bigP = BigInteger.valueOf(p);
        while(N.mod(bigP) == BigInteger.ZERO) {
            //System.out.println("N: "+ N + ", N%p: " + N % p);
            N = N.divide(bigP);
            i++;
        }
        int res = i%2;
        if (res == 0) {
            return false;
        } else {
            return true;
        }
    }

    public static int hasEvenNumberOfPrime(int N, int p) {
        int i = 0;
        while(N % p == 0) {
            //System.out.println("N: "+ N + ", N%p: " + N % p);
            N = N/p;
            i++;
        }
        return i%2;
    }

    private static int gcdThing(int a, int b) {
        BigInteger b1 = new BigInteger(""+a); // there's a better way to do this. I forget.
        BigInteger b2 = new BigInteger(""+b);
        BigInteger gcd = b1.gcd(b2);
        return gcd.intValue();
    }

    private static boolean[] gaussWrapper(boolean[][] inputMatrix, int test) {
        int bitlength = inputMatrix.length;
        BitSet[] bitarray = gaussGF2(inputMatrix);
        int[][] convertedMatrix = convertMatrix(bitarray, bitlength);
        boolean[] dependencies = getDependentRows(convertedMatrix, test);

        return dependencies;
    }

    private static BitSet[] gaussGF2(boolean[][] inputMatrix) {
        BitSet[] bitArray = new BitSet[inputMatrix[0].length];
        for(int i = 0; i < inputMatrix[0].length; i++) {
            bitArray[i] = new BitSet(inputMatrix.length);
            for(int j = 0; j < inputMatrix.length; j++) {
                if(inputMatrix[j][i]) {
                    bitArray[i].set(j);
                }
            }
        }
        markedRows = new boolean[inputMatrix.length];
        for(int col = 0; col < bitArray.length; col++) {
            int nextSetBit = bitArray[col].nextSetBit(0);
            if(nextSetBit == -1) continue;
            markedRows[nextSetBit] = true;
            for(int c = 0; c < bitArray.length; c++) {
                if(c == col) continue;
                if(bitArray[c].get(nextSetBit)) bitArray[c].xor(bitArray[col]);
            }
        }
        return bitArray;
    }

    private static boolean[] getDependentRows(int[][] matrix, int test) {
        int count = 0;
        boolean[] dependent = new boolean[matrix.length];
        for(int i = 0; i < matrix.length; i++) {
            if(markedRows[i]) continue;
            if(count != test) {
                count++;
                continue;
            }
            for(int j = 0; j < matrix[0].length; j++) {
                if(matrix[i][j] == 0) continue;
                for(int row = 0; row < matrix.length; row++)
                    if(matrix[row][j] == 1) dependent[row] = true;
            }
            break;
        }
        return dependent;
    }
    private static int[][] convertMatrix(BitSet[] bitArray, int bitlength) {
        int[][] matrix = new int[bitlength][bitArray.length];
        for(int i = 0; i < bitlength; i++) {
            for(int j = 0; j < bitArray.length; j++) {
                if(bitArray[j].get(i))
                    matrix[i][j] = 1;
                else
                    matrix[i][j] = 0;
            }
        }
        return matrix;
    }

}