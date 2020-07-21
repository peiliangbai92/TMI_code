// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>



using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
mat soft_full(mat L, double lambda){
	int height=L.n_rows;
	int width=L.n_cols;
 	for(int y = 0; y < height; y++){
 		for(int x=0; x < width; x++) {
            if (L(y, x) > lambda)
            	L(y, x) = L(y, x) - lambda;
            else if (L(y, x) <  - lambda)
                L(y, x) = L(y, x) + lambda;
            else
                L(y, x) = 0.0;
        }

 	} 
 	return L;
}

// [[Rcpp::export]]
List var_break_fit_block_cpp(NumericMatrix data, double lambda, double lambda2, int p, int max_iteration, double tol , NumericMatrix initial_phi, NumericVector blocks, NumericVector cv_index ){

    int k = data.ncol(); int T = data.nrow(); int n_new = blocks.size() - 1;
    int n = T - p;

    // Rprintf("%i", k);
    // Rprintf("%i", T);
    // Rprintf("%i", n_new);
    // Rprintf("%i", n);

    mat data_m(data.begin(), T, k);
    List Y_b(n_new);
    List y_b(n_new);
    //Rcout  << size(data_m)[1];

    for(int i = 0; i < n_new; i++) {
    	y_b[i] = data_m( span(blocks[i],  blocks[i+1]-1 ), span::all) ;
    }
    y_b[0] = data_m(  span(p, (blocks[1]-1) ), span::all );


    mat Y( k*p, n);
    for(int i = (p-1); i < (T-1); i++) {
    	for(int j = 1; j <= p; j++){
    		Y.submat( (j-1)*k,(i-p+1), j*k-1 ,(i-p+1) ) = data_m.submat(i-j+1, 0, i-j+1,  k-1 ).t();
    	}   	

    }
    
    Y_b[0] = Y.submat( 0 ,0, k*p-1,blocks[1]-p-1  );
    for(int i = 1; i < n_new; i++) {
    	Y_b[i] = Y.submat( 0, blocks[i]-p ,k*p-1 , blocks[i+1]-p-1  ) ;
    }
    
    //NumericMatrix temp = y_b[0]; 
    //Rcout << temp.nrow();

    int cv_l = cv_index.size();
    if( cv_l >0){
    	for(int t_1 = 0; t_1 < cv_l; t_1 ++){
    		mat yb_temp = y_b[cv_index[t_1]-1]; 
    		mat Yb_temp = Y_b[cv_index[t_1]-1]; 
    		int tt =  yb_temp.n_rows;
    		//Rcout << tt1;
            y_b[cv_index[t_1]-1] = yb_temp(span(0,tt-2), span::all );
            Y_b[cv_index[t_1]-1] = Yb_temp( span::all ,span(0,tt-2));

    	}
    } 

    List C(n_new);
    for(int j = 0; j < n_new; j++){
    	mat yb_temp = y_b[j];
    	mat Yb_temp = Y_b[j]; 
    	C[j] = Yb_temp * yb_temp;
    }  
    
    mat C_sum (k*p*n_new, k, fill::zeros);
    //Rcout << C_sum(1,1);
    mat C_temp = C[0];
    C_sum(span(0, k*p - 1) , span::all) = C_temp;
    for(int i = 2; i <= n_new ; i++){
    	mat C_temp = C[i-1];
    	//Rcout <<  size(C_temp);
    	C_sum(  span((i-1)*k*p, (i*k*p) -1  ), span::all ) = C_sum( span( (i-2)*k*p,(i-1)*k*p -1 ), span::all) + C_temp; 
    }
    mat C_sum_new (k*p*n_new,k, fill::zeros);
    C_sum_new( span(   0, k*p - 1 ), span::all) =  C_sum(   span(  (n_new-1)*k*p,  n_new*k*p - 1  ), span::all );
    for(int i = 2; i <= n_new ; i++){
    	C_sum_new(  span((i-1)*k*p, (i*k*p) -1  ), span::all )  = C_sum(   span(  (n_new-1)*k*p,  n_new*k*p - 1  ), span::all ) - C_sum( span ((i-2)*k*p, (i-1)*k*p -1 ), span::all);
    }


    List D(n_new);
    for(int j = 0; j < n_new; j++){
    	mat Yb_temp = Y_b[j]; 
    	D[j] = Yb_temp * Yb_temp.t();
    }

    mat D_sum (k*p*n_new,k*p, fill::zeros);
    mat D_temp = D[0];
    D_sum(span(0, k*p - 1) , span::all) = D_temp;
    for(int i = 2; i <= n_new ; i++){
    	mat D_temp = D[i-1];
    	//Rcout <<  size(C_temp);
    	D_sum(  span((i-1)*k*p, (i*k*p) -1  ), span::all ) = D_sum( span( (i-2)*k*p,(i-1)*k*p -1 ), span::all) + D_temp; 
    }
    mat D_sum_new (k*p*n_new,k*p, fill::zeros);
    D_sum_new( span(  0, k*p - 1 ), span::all) =  D_sum(   span(  (n_new-1)*k*p,  n_new*k*p - 1  ), span::all );
    for(int i = 2; i <= n_new ; i++){
    	D_sum_new(  span((i-1)*k*p, (i*k*p) -1  ), span::all )  = D_sum(   span(  (n_new-1)*k*p,  n_new*k*p - 1  ), span::all ) - D_sum( span ((i-2)*k*p, (i-1)*k*p -1 ), span::all);
    }

    mat D_sum_new_inv (k*p*n_new,k*p, fill::zeros);
    // for(int i = 1; i <= n_new ; i++){
    // 	mat noise(k*p,k*p, fill::eye);
    // 	D_sum_new_inv(  span((i-1)*k*p, (i*k*p) -1  ), span::all )  = inv( D_sum_new( span((i-1)*k*p, (i*k*p) -1  ), span::all  ) +  pow(10, -6)*noise );
    // }

    for(int i = 1; i <= n_new ; i++){
        vec eigval =  eig_sym(D_sum_new ( span((i-1)*k*p, (i*k*p)-1), span::all)  );
        double add_pd =0.0;
        double min_eigen = eigval(0);
        //Rcout << eigval << '\n';
        //Rcout << min_eigen << '\n';
        if(min_eigen <= 0){
            Rprintf("Invertiable! adding noise!");
            add_pd = (10)*fabs( min_eigen);
        }
        mat noise(k*p,k*p, fill::eye);
        D_sum_new_inv(  span((i-1)*k*p, (i*k*p) -1  ), span::all )  = inv( D_sum_new( span((i-1)*k*p, (i*k*p) -1  ), span::all  ) +  add_pd*noise );
    }

    // mat D_new (n_new, p*k, fill::zeros);
    // for(int i = 0; i < n_new ; i++){
    // 	mat Yb_temp = Y_b[i]; 
    // 	mat D_new_temp  = diagvec(Yb_temp * Yb_temp.t() ); 
    // 	D_new(i , span::all) =	reshape(D_new_temp, 1, p*k);
    // }


    mat phi_hat(initial_phi.begin(), k, k*p*n_new);
    mat phi_new (k, k*p*n_new, fill::zeros);

    int flag = 0;
    int l = 2;
    //Rcout <<  floor(0.5 * max_iteration); 
    //Rcout <<  (l == floor(0.5 * max_iteration)); 
    
    while( l < max_iteration){
    	if(  l == floor(0.5 * max_iteration)  ){
    		tol = (2)*tol;
    	}
    	if( l  == floor( 0.75 * max_iteration )) {
    		tol = (4/2)*tol;
    	}
        l = l+1; 
        mat phi_compare = phi_hat;

        for(int i = 1; i <= n_new ; i++){
        	List E(n_new);
        	for(int j = 1; j <= n_new ; j++){
        		
        		E[j-1] = D_sum_new( span(  ( std::max( j,i )-1)*k*p  , std::max(j,i )*k*p  -1 ), span::all)  * phi_hat(span::all , span( (j-1)*k*p, j*k*p -1  )).t() ;    			
        	}

        	
        	

        	mat E_new = E[0];
    		for ( int g = 1; g < n_new; g++ ) {
    			mat E_temp = E[g];
    			E_new = E_new + E_temp;
    		}


    		

    		E_new =  E_new - D_sum_new(  span( (i-1)*k*p,  i*k*p -1  ), span::all)  * phi_hat( span::all, span( (i-1)*k*p, i*k*p -1) ).t();

    		mat S =  C_sum_new(  span(  (i-1)*k*p , i*k*p -1 ), span::all )  - E_new;




    		

    		S = soft_full(S, lambda);
            // Rcout << S;

    		// if(i == 1){
      //   		mat test_temp = S;
      //   		Rcout << test_temp;
      //   	}

    		mat phi_temp = D_sum_new_inv( span ((i-1)*k*p , i*k*p-1), span::all )  *  S;
    		
    		// mat prod (k*p,k*p, fill::zeros);
    		// for(int w = 0; w < k*p; w++){
    		// 	prod(w , w ) = 1/sum(D_new(  span(i-1, n_new-1), w )) ;  
    		// }
    		phi_temp = phi_temp.t();

            // Rprintf("i= ");
            // Rcout <<  i;
            // Rprintf("\n ");
            // Rcout << phi_temp;
            // if(i == 31){
            //     mat test_temp = S;
            //     Rcout << test_temp;
            //     Rprintf("\n ");
            //     Rcout << phi_temp;
            // }

    		phi_hat( span::all ,   span( (i-1)*k*p, i*k*p -1)  ) = phi_temp;
            phi_new( span::all ,   span( (i-1)*k*p, i*k*p -1)  ) = phi_temp;
                                                                  

        }

        //Rcout<< phi_new( span::all ,   span( k*p, 2*k*p -1)  ) ;


        mat phi_temp_soft(k,k*p*n_new, fill::zeros);
        phi_temp_soft( span::all , span(0, k*p-1) ) = soft_full(phi_hat( span::all ,  span(0, k*p-1) ),lambda2);
        mat temp_1 = phi_hat( span::all ,  span(0, k*p-1) );
        for(int z_1 = 2; z_1 <= n_new; z_1 ++){
            mat temp_2 = temp_1 + phi_hat( span::all ,  span((z_1-1)*k*p, z_1*k*p-1) );
            mat temp_1_soft = soft_full(temp_1,lambda2);
            mat temp_2_soft = soft_full(temp_2,lambda2);
            phi_temp_soft(   span::all ,  span( (z_1-1)*k*p , z_1*k*p  -1 ) )= temp_2_soft - temp_1_soft;
            temp_1 = temp_2;
        }
        phi_new  = phi_temp_soft;



        mat abs_temp = abs(phi_new - phi_compare);
        //Rcout << max_temp;
        double max_temp = abs_temp.max();
        //Rcout << max_temp;
        //Rcout << '\n';
        if ( max_temp < tol) {
            break;
        } 
        if (  max_temp > pow(10, 5)) {
            Rprintf("NOT CONVERGED");
            flag = 1;
            break;
        }
        if ( max_temp > tol ) {
        	phi_hat = phi_new;    
            Rprintf( "%f \n", max_temp);
        }                
	}

	

    //return List::create(Named("phi.hat")= phi_hat,Named("flag")= flag, Named("D_sum_new_inv") = D_sum_new_inv, Named("D_sum_new") = D_sum_new);
    return List::create(Named("phi.hat")= phi_hat,Named("flag")= flag);

	
}



// [[Rcpp::export]]
List lambda_warm_up(NumericMatrix data, int p, NumericVector blocks, NumericVector cv_index ){

    int k = data.ncol(); int T = data.nrow(); int n_new = blocks.size() - 1;
    int n = T - p;


    mat data_m(data.begin(), T, k);
    List Y_b(n_new);
    List y_b(n_new);

    for(int i = 0; i < n_new; i++) {
        y_b[i] = data_m( span(blocks[i],  blocks[i+1]-1 ), span::all) ;
    }
    y_b[0] = data_m(  span(p, (blocks[1]-1) ), span::all );


    mat Y( k*p, n);
    for(int i = (p-1); i < (T-1); i++) {
        for(int j = 1; j <= p; j++){
            Y.submat( (j-1)*k,(i-p+1), j*k-1 ,(i-p+1) ) = data_m.submat(i-j+1, 0, i-j+1,  k-1 ).t();
        }       

    }
    
    Y_b[0] = Y.submat( 0 ,0, k*p-1,blocks[1]-p-1  );
    for(int i = 1; i < n_new; i++) {
        Y_b[i] = Y.submat( 0, blocks[i]-p ,k*p-1 , blocks[i+1]-p-1  ) ;
    }
    

    int cv_l = cv_index.size();
    if( cv_l >0){
        for(int t_1 = 0; t_1 < cv_l; t_1 ++){
            mat yb_temp = y_b[cv_index[t_1]-1]; 
            mat Yb_temp = Y_b[cv_index[t_1]-1]; 
            int tt =  yb_temp.n_rows;
            //Rcout << tt1;
            y_b[cv_index[t_1]-1] = yb_temp(span(0,tt-2), span::all );
            Y_b[cv_index[t_1]-1] = Yb_temp( span::all ,span(0,tt-2));

        }
    } 

    List C(n_new);
    for(int j = 0; j < n_new; j++){
        mat yb_temp = y_b[j];
        mat Yb_temp = Y_b[j]; 
        C[j] = Yb_temp * yb_temp;
    }  
    
    mat C_sum (k*p*n_new, k, fill::zeros);
    mat C_temp = C[0];
    C_sum(span(0, k*p - 1) , span::all) = C_temp;
    for(int i = 2; i <= n_new ; i++){
        mat C_temp = C[i-1];
        C_sum(  span((i-1)*k*p, (i*k*p) -1  ), span::all ) = C_sum( span( (i-2)*k*p,(i-1)*k*p -1 ), span::all) + C_temp; 
    }
    mat C_sum_new (k*p*n_new,k, fill::zeros);
    C_sum_new( span(   0, k*p - 1 ), span::all) =  C_sum(   span(  (n_new-1)*k*p,  n_new*k*p - 1  ), span::all );
    for(int i = 2; i <= n_new ; i++){
        C_sum_new(  span((i-1)*k*p, (i*k*p) -1  ), span::all )  = C_sum(   span(  (n_new-1)*k*p,  n_new*k*p - 1  ), span::all ) - C_sum( span ((i-2)*k*p, (i-1)*k*p -1 ), span::all);
    }

    double lambda_max = 0;
    for(int i = 1; i <= n_new; i++){
        double lambda_temp = abs(C_sum_new(  span((i-1)*k*p, (i*k*p) -1  ), span::all )).max(); 
        // Rcout << lambda_temp;
        // Rcout << '\n';
        // Rcout << lambda_max;
        // Rcout << '\n';
        // Rcout << std::max(lambda_max,lambda_temp);
        // Rcout << '\n';
        lambda_max = std::max(lambda_max,lambda_temp);

    }


    return List::create(Named("lambda_1_max")= lambda_max);

    
}

