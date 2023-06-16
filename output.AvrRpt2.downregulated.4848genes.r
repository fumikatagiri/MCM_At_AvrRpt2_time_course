load('./data/AvrRpt2_negative_genes.Rdata')
nega.df = data.frame(AvrRpt2.nega.genes = AvrRpt2_mock_negative_genes)
write.csv(nega.df, file = './data/AvrRpt2.downregulated.genes.4848.csv')
