using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ICD.Beagle
{
    class Genome
    {
        public List<Gene> Genes;
        public double Fitness;

        public Genome()
        {
            Genes = new List<Gene>();
            Fitness = double.NaN;
        }

        public Genome(List<double> geneMinValues, List<double> geneMaxValues, List<double> geneValues = null)
        {
            Genes = new List<Gene>();

            if (geneValues == null)
                for (int i = 0; i < geneMinValues.Count; i++)
                    Genes.Add(new Gene(geneMinValues[i], geneMaxValues[i]));
            else
                for (int i = 0; i < geneMinValues.Count; i++)
                    Genes.Add(new Gene(geneMinValues[i], geneMaxValues[i], geneValues[i]));
        }


        public double DistanceTo(Genome targetGenome)
        {
            double distanceSquared = 0.0;
            for (int i = 0; i < Genes.Count; i++)
                distanceSquared += Math.Pow(Genes[i].Value - targetGenome.Genes[i].Value, 2.0);

            return Math.Sqrt(distanceSquared);
        }

        public void MakeRandom()
        {
            foreach (Gene gene in Genes) gene.MakeRandom();
        }

        public List<double> GetGeneValues()
        {
            List<double> genes = new List<double>();
            foreach (Gene gene in this.Genes)
                genes.Add(gene.Value);
            return genes;
        }
    }
}
