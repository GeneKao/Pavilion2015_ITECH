using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ICD.Beagle
{
    class Population
    {
        public List<Genome> Genomes;

        public Population()
        {
            Genomes = new List<Genome>();
        }

        public Population(List<double> geneMinValues, List<double> geneMaxValues)
        {
            Genomes = new List<Genome>();
            for (int i = 0; i < geneMinValues.Count; i++)
                Genomes.Add(new Genome(geneMinValues, geneMaxValues));
        }

        public void MakeRandom()
        {
            foreach (Genome genome in Genomes)
                genome.MakeRandom();
        }
    }
}
