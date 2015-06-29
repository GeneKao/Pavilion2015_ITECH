using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ICD;

namespace ICD.Beagle
{
    class Beagle
    {
        public int PopulationSize = 0;
        public List<Genome> Genomes = null;
        public double MaxCouplingDistance = double.NaN;
        public double MutationRate = double.NaN;
        public double MutationAmount = double.NaN;
        public ComputeFitnessDelegate computeFitness = null;

        private List<Genome> parents = null;      

        public Beagle()
        {
        }


        public Beagle(int populationSize, List<double> geneMinValues, List<double> geneMaxValues, ComputeFitnessDelegate computeFitness = null)
        {
            this.PopulationSize = populationSize;

            Genomes = new List<Genome>();

            for (int i = 0; i < populationSize; i++)
                Genomes.Add(new Genome(geneMinValues, geneMaxValues));

            parents = new List<Genome>();

            this.computeFitness = computeFitness;
            if (computeFitness != null)
                foreach (Genome genome in Genomes)
                    computeFitness(genome);
        }


        public void Evolve()
        {
            SelectParents();
            Reproduce();
            SelectFittests();
        }


        public void SelectParents()
        {
            parents.Clear();

            for (int i = 0; i < PopulationSize; i++)
            {
                List<int> potentialMateIndices = new List<int>();
                for (int j = 0; j < PopulationSize; j++)
                    if (i != j && Genomes[i].DistanceTo(Genomes[j]) < MaxCouplingDistance) potentialMateIndices.Add(j);

                if (potentialMateIndices.Count > 0)
                {
                    parents.Add(Genomes[i]);
                    parents.Add(Genomes[potentialMateIndices[Utils.GetRandomInteger(potentialMateIndices.Count - 1)]]);
                }             
            }
        }


        public void Reproduce()
        {
            int N = (int)(parents.Count * 0.5);

            for (int i = 0; i < N; i++)
            {
                makeChild(parents[i * 2], parents[i * 2 + 1]);
            }
        }


        private void makeChild(Genome father, Genome mother)
        {
            Genome child = new Genome();
            for (int i = 0; i < father.Genes.Count; i++)
            {
                Gene fatherGene = father.Genes[i];
                Gene motherGene = mother.Genes[i];
                Gene childGene = new Gene(fatherGene.Min, fatherGene.Max, 0.5 * (fatherGene.Value + motherGene.Value));

                if (Utils.GetRandomDouble() < MutationRate)
                {
                    childGene.Value = childGene.Value + (childGene.Max - childGene.Min) * Utils.GetRandomDouble(-MutationAmount,  MutationAmount);
                    childGene.MakeValueValid();
                }

                child.Genes.Add(childGene);
            }

            //this.computeFitness(child);
            Genomes.Add(child);
        }


        public void SelectFittests()
        {
            Genomes.Sort((A, B) => (A.Fitness < B.Fitness ? -1 : (A.Fitness > B.Fitness ? 1 : 0)));
            Genomes.RemoveRange(PopulationSize, Genomes.Count - PopulationSize);
        }
    }
}
