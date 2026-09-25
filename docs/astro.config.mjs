import starlight from '@astrojs/starlight';
import { defineConfig } from 'astro/config';
import starlightPydocs, { pydocsSidebarGroup } from 'starlight-pydocs';

const REPO = 'https://github.com/vagkaratzas/mgnifam';

export default defineConfig({
  site: 'https://vagkaratzas.github.io',
  base: '/mgnifam',
  integrations: [
    starlight({
      title: 'mgnifam',
      description: 'Iterative HMM-based protein family generation over very large sequence databases.',
      favicon: '/favicon.png',
      social: [{ icon: 'github', label: 'GitHub', href: REPO }],
      editLink: { baseUrl: `${REPO}/edit/main/docs/` },
      plugins: [
        starlightPydocs({
          packages: [
            {
              name: 'mgnifam',
              search: ['../src'],
              sourceLink: { host: 'github', repo: 'vagkaratzas/mgnifam', ref: 'main', root: '..' },
            },
          ],
          inventories: ['python'],
        }),
      ],
      sidebar: [
        { label: 'Home', link: '/' },
        {
          label: 'Guides',
          items: [
            { label: 'Installation', link: '/guides/installation/' },
            { label: 'Generating families', link: '/guides/generate-families/' },
            { label: 'Updating families', link: '/guides/update-families/' },
          ],
        },
        {
          label: 'Reference',
          items: [
            { label: 'Outputs', link: '/reference/outputs/' },
            { label: 'Reproducibility', link: '/reference/reproducibility/' },
            { label: 'Development', link: '/reference/development/' },
            { label: 'Changelog', link: `${REPO}/blob/main/CHANGELOG.md` },
          ],
        },
        {
          label: 'Python API',
          items: [{ label: 'Stability', link: '/reference/python-api/' }, pydocsSidebarGroup],
        },
      ],
    }),
  ],
});
