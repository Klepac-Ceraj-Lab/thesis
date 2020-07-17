library(rotl)
taxa<-tnrs_match_names(names= c("Danaus plexippus",
                                "Sylvia atricapilla",
                                "Acrocephalus arundinaceus",
                                "Anas platyrhynchos",
                                "Oncorhynchus nerka",
                                "Hirundo pyrrhonota",
                                "Progne subis"))

tree <- tol_induced_subtree(ott_ids = ott_id(taxa))

plot(tree, cex = .8, label.offset = .1, no.margin = TRUE)