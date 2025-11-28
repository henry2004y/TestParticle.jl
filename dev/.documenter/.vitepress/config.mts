import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";
import path from 'path'

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: '/henry2004y.github.io/TestParticle.jl/dev/',// TODO: replace this in makedocs!
}

const navTemp = {
  nav: [
{ text: 'Home', link: '/index' },
{ text: 'Tutorial', link: '/tutorial' },
{ text: 'Examples', collapsed: false, items: [
{ text: 'Overview', link: '/generated/index' },
{ text: 'Basics', collapsed: false, items: [
{ text: 'generated/basics/demo_energy_conservation', link: '/generated/basics/demo_energy_conservation' },
{ text: 'generated/basics/demo_boris', link: '/generated/basics/demo_boris' },
{ text: 'generated/basics/demo_Buniform_Ezero', link: '/generated/basics/demo_Buniform_Ezero' },
{ text: 'generated/basics/demo_dimensionless', link: '/generated/basics/demo_dimensionless' },
{ text: 'generated/basics/demo_dimensionless_periodic', link: '/generated/basics/demo_dimensionless_periodic' },
{ text: 'generated/basics/demo_dimensionless_dimensional', link: '/generated/basics/demo_dimensionless_dimensional' },
{ text: 'generated/basics/demo_electron_proton', link: '/generated/basics/demo_electron_proton' },
{ text: 'generated/basics/demo_multiple', link: '/generated/basics/demo_multiple' },
{ text: 'generated/basics/demo_ExB_drift', link: '/generated/basics/demo_ExB_drift' },
{ text: 'generated/basics/demo_gravity_drift', link: '/generated/basics/demo_gravity_drift' },
{ text: 'generated/basics/demo_gradient_B', link: '/generated/basics/demo_gradient_B' },
{ text: 'generated/basics/demo_curvature_B', link: '/generated/basics/demo_curvature_B' },
{ text: 'generated/basics/demo_FLR', link: '/generated/basics/demo_FLR' },
{ text: 'generated/basics/demo_polarization_drift', link: '/generated/basics/demo_polarization_drift' },
{ text: '–-', link: '/generated/basics/demo_array' }]
 },
{ text: 'Advanced', collapsed: false, items: [
{ text: 'generated/advanced/demo_boris_outofdomain', link: '/generated/advanced/demo_boris_outofdomain' },
{ text: 'generated/advanced/demo_cosmicray', link: '/generated/advanced/demo_cosmicray' },
{ text: 'generated/advanced/demo_ensemble', link: '/generated/advanced/demo_ensemble' },
{ text: 'generated/advanced/demo_flux', link: '/generated/advanced/demo_flux' },
{ text: 'generated/advanced/demo_savingcallback', link: '/generated/advanced/demo_savingcallback' },
{ text: 'generated/advanced/demo_output_func', link: '/generated/advanced/demo_output_func' },
{ text: 'generated/advanced/demo_currentsheet', link: '/generated/advanced/demo_currentsheet' },
{ text: 'generated/advanced/demo_magneticmirror', link: '/generated/advanced/demo_magneticmirror' },
{ text: 'generated/advanced/demo_magneticbottle', link: '/generated/advanced/demo_magneticbottle' },
{ text: 'generated/advanced/demo_proton_dipole', link: '/generated/advanced/demo_proton_dipole' },
{ text: 'generated/advanced/demo_analytic_magnetosphere', link: '/generated/advanced/demo_analytic_magnetosphere' },
{ text: 'generated/advanced/demo_shock', link: '/generated/advanced/demo_shock' },
{ text: 'generated/advanced/demo_fermi_foreshock', link: '/generated/advanced/demo_fermi_foreshock' },
{ text: 'generated/advanced/demo_spherical', link: '/generated/advanced/demo_spherical' },
{ text: 'generated/advanced/demo_tokamak_coil', link: '/generated/advanced/demo_tokamak_coil' },
{ text: 'generated/advanced/demo_tokamak_profile', link: '/generated/advanced/demo_tokamak_profile' },
{ text: 'generated/advanced/demo_gc', link: '/generated/advanced/demo_gc' },
{ text: '–-', link: '/generated/advanced/demo_batsrus_3dstructured' },
{ text: 'generated/advanced/demo_radiation', link: '/generated/advanced/demo_radiation' },
{ text: '–-', link: '/generated/advanced/demo_gpu' }]
 }]
 },
{ text: 'API', link: '/api' },
{ text: 'Plot Functions', link: '/plotfunctions' }
]
,
}

const nav = [
  ...navTemp.nav,
  {
    component: 'VersionPicker'
  }
]

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: '/henry2004y.github.io/TestParticle.jl/dev/',// TODO: replace this in makedocs!
  title: 'TestParticle.jl',
  description: 'Test particle tracing in fields',
  lastUpdated: true,
  cleanUrls: true,
  ignoreDeadLinks: true,
  outDir: '../1', // This is required for MarkdownVitepress to work correctly...
  head: [
    
    ['script', {src: `${getBaseRepository(baseTemp.base)}versions.js`}],
    // ['script', {src: '/versions.js'], for custom domains, I guess if deploy_url is available.
    ['script', {src: `${baseTemp.base}siteinfo.js`}]
  ],
  
  vite: {
    define: {
      __DEPLOY_ABSPATH__: JSON.stringify('/henry2004y.github.io/TestParticle.jl'),
    },
    resolve: {
      alias: {
        '@': path.resolve(__dirname, '../components')
      }
    },
    optimizeDeps: {
      exclude: [ 
        '@nolebase/vitepress-plugin-enhanced-readabilities/client',
        'vitepress',
        '@nolebase/ui',
      ], 
    }, 
    ssr: { 
      noExternal: [ 
        // If there are other packages that need to be processed by Vite, you can add them here.
        '@nolebase/vitepress-plugin-enhanced-readabilities',
        '@nolebase/ui',
      ], 
    },
  },
  markdown: {
    math: true,
    config(md) {
      md.use(tabsMarkdownPlugin),
      md.use(mathjax3),
      md.use(footnote)
    },
    theme: {
      light: "github-light",
      dark: "github-dark"}
  },
  themeConfig: {
    outline: 'deep',
    
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar: [
{ text: 'Home', link: '/index' },
{ text: 'Tutorial', link: '/tutorial' },
{ text: 'Examples', collapsed: false, items: [
{ text: 'Overview', link: '/generated/index' },
{ text: 'Basics', collapsed: false, items: [
{ text: 'generated/basics/demo_energy_conservation', link: '/generated/basics/demo_energy_conservation' },
{ text: 'generated/basics/demo_boris', link: '/generated/basics/demo_boris' },
{ text: 'generated/basics/demo_Buniform_Ezero', link: '/generated/basics/demo_Buniform_Ezero' },
{ text: 'generated/basics/demo_dimensionless', link: '/generated/basics/demo_dimensionless' },
{ text: 'generated/basics/demo_dimensionless_periodic', link: '/generated/basics/demo_dimensionless_periodic' },
{ text: 'generated/basics/demo_dimensionless_dimensional', link: '/generated/basics/demo_dimensionless_dimensional' },
{ text: 'generated/basics/demo_electron_proton', link: '/generated/basics/demo_electron_proton' },
{ text: 'generated/basics/demo_multiple', link: '/generated/basics/demo_multiple' },
{ text: 'generated/basics/demo_ExB_drift', link: '/generated/basics/demo_ExB_drift' },
{ text: 'generated/basics/demo_gravity_drift', link: '/generated/basics/demo_gravity_drift' },
{ text: 'generated/basics/demo_gradient_B', link: '/generated/basics/demo_gradient_B' },
{ text: 'generated/basics/demo_curvature_B', link: '/generated/basics/demo_curvature_B' },
{ text: 'generated/basics/demo_FLR', link: '/generated/basics/demo_FLR' },
{ text: 'generated/basics/demo_polarization_drift', link: '/generated/basics/demo_polarization_drift' },
{ text: '–-', link: '/generated/basics/demo_array' }]
 },
{ text: 'Advanced', collapsed: false, items: [
{ text: 'generated/advanced/demo_boris_outofdomain', link: '/generated/advanced/demo_boris_outofdomain' },
{ text: 'generated/advanced/demo_cosmicray', link: '/generated/advanced/demo_cosmicray' },
{ text: 'generated/advanced/demo_ensemble', link: '/generated/advanced/demo_ensemble' },
{ text: 'generated/advanced/demo_flux', link: '/generated/advanced/demo_flux' },
{ text: 'generated/advanced/demo_savingcallback', link: '/generated/advanced/demo_savingcallback' },
{ text: 'generated/advanced/demo_output_func', link: '/generated/advanced/demo_output_func' },
{ text: 'generated/advanced/demo_currentsheet', link: '/generated/advanced/demo_currentsheet' },
{ text: 'generated/advanced/demo_magneticmirror', link: '/generated/advanced/demo_magneticmirror' },
{ text: 'generated/advanced/demo_magneticbottle', link: '/generated/advanced/demo_magneticbottle' },
{ text: 'generated/advanced/demo_proton_dipole', link: '/generated/advanced/demo_proton_dipole' },
{ text: 'generated/advanced/demo_analytic_magnetosphere', link: '/generated/advanced/demo_analytic_magnetosphere' },
{ text: 'generated/advanced/demo_shock', link: '/generated/advanced/demo_shock' },
{ text: 'generated/advanced/demo_fermi_foreshock', link: '/generated/advanced/demo_fermi_foreshock' },
{ text: 'generated/advanced/demo_spherical', link: '/generated/advanced/demo_spherical' },
{ text: 'generated/advanced/demo_tokamak_coil', link: '/generated/advanced/demo_tokamak_coil' },
{ text: 'generated/advanced/demo_tokamak_profile', link: '/generated/advanced/demo_tokamak_profile' },
{ text: 'generated/advanced/demo_gc', link: '/generated/advanced/demo_gc' },
{ text: '–-', link: '/generated/advanced/demo_batsrus_3dstructured' },
{ text: 'generated/advanced/demo_radiation', link: '/generated/advanced/demo_radiation' },
{ text: '–-', link: '/generated/advanced/demo_gpu' }]
 }]
 },
{ text: 'API', link: '/api' },
{ text: 'Plot Functions', link: '/plotfunctions' }
]
,
    editLink: { pattern: "https://https://github.com/henry2004y/TestParticle.jl/edit/master/docs/src/:path" },
    socialLinks: [
      { icon: 'github', link: 'https://github.com/henry2004y/TestParticle.jl' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/dev/" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`
    }
  }
})
