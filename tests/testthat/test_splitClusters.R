context('splitClusters')

ttree <- "(((t1:0.2,t8:0.15):0.05,(t5:0.1,t4:0.1):0.2):0.2,(t6:0.2,(t3:0.1,(t2:0.07,t7:0.07):0.05):0.15):0.1);"
phy <- ape::read.tree(text = ttree)
phy$tip.label <- c('genome1;gene1','genome1;gene2','genome2;gene1','genome3;gene1',
                   'genome1;gene3','genome2;gene2','genome3;gene2','genome3;gene3')



l <- ape::subtrees(phy)
x <- pewit:::whichNodesHaveRecentParalogues(subtrees = l)


test_that('whichNodesHaveRecentParalogues work', {
  expect_is(x, 'integer')
  expect_length(x, 2L)
  expect_equal(x, c(11L, 15L))
})


d <- as.dist(ape::cophenetic.phylo(phy))
tree.trim <- pewit:::rmvRecentParalogues(phy,w = x,d = d)
test_that('rmvRecentParalogues works',{
  expect_is(tree.trim, 'list')
  expect_named(tree.trim, c('tree', 'paralogs'))
  expect_is(tree.trim$tree, 'phylo')
  expect_is(tree.trim$paralogs, 'list')
  expect_equal(tree.trim$tree$Nnode, 5)
  expect_length(tree.trim$tree$tip.label, 6)
  expect_length(tree.trim$paralogs, 2)

})

l <- ape::subtrees(tree.trim$tree)
x <- pewit:::whichSubtreesAreTrueOrthologues(subtrees = l)
test_that('whichSubtreesAreTrueOrthologues works',{
  expect_is(x, 'logical')
  expect_length(x, 5L)
  expect_equal(x, c(FALSE, TRUE, TRUE, TRUE, TRUE))
})

y <- pewit:::maxTrueOGs(tru.or = x, l = l)
test_that('maxTrueOGs works',{
  expect_is(y, 'logical')
  expect_length(y, 5L)
  expect_equal(y, c(FALSE, TRUE, FALSE, TRUE, FALSE))
})


